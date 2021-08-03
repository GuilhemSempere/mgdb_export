/*******************************************************************************
 * MGDB Export - Mongo Genotype DataBase, export handlers
 * Copyright (C) 2016 - 2019, <CIRAD> <IRD>
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License, version 3 as published by
 * the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
 * details.
 *
 * See <http://www.gnu.org/licenses/agpl.html> for details about GNU General
 * Public License V3.
 *******************************************************************************/
package fr.cirad.mgdb.exporting.markeroriented;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.OutputStream;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.bson.Document;
import org.springframework.data.mongodb.core.MongoTemplate;

import com.mongodb.client.MongoCollection;

import fr.cirad.mgdb.exporting.AbstractExportWritingThread;
import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.exporting.tools.ExportManager;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingSample;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
import fr.cirad.mgdb.model.mongo.subtypes.SampleGenotype;
import fr.cirad.tools.AlphaNumericComparator;
import fr.cirad.tools.Helper;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mongo.MongoTemplateManager;
import htsjdk.variant.variantcontext.VariantContext.Type;

/**
 * The Class GFFExportHandler.
 */
public class GFFExportHandler extends AbstractMarkerOrientedExportHandler {

    /**
     * The Constant LOG.
     */
    private static final Logger LOG = Logger.getLogger(GFFExportHandler.class);

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatName()
     */
    @Override
    public String getExportFormatName() {
        return "GFF3";
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
     */
    @Override
    public String getExportFormatDescription() {
        return "Exports data in GFF3 Format based on Sequence Ontology. See <a target='_blank' href='https://m.ensembl.org/info/website/upload/gff3.html'>https://m.ensembl.org/info/website/upload/gff3.html</a> and <a target='_blank' href='https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md'>https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md</a>";
    }
    
	@Override
	public String getExportArchiveExtension() {
		return "zip";
	}

    @Override
    public void exportData(OutputStream outputStream, String sModule, Collection<String> individuals1, Collection<String> individuals2, ProgressIndicator progress, String tmpVarCollName, Document varQuery, long markerCount, Map<String, String> markerSynonyms, HashMap<String, Float> annotationFieldThresholds, HashMap<String, Float> annotationFieldThresholds2, List<GenotypingSample> samplesToExport, Collection<String> individualMetadataFieldsToExport, Map<String, InputStream> readyToExportFiles) throws Exception { 
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        ZipOutputStream zos = IExportHandler.createArchiveOutputStream(outputStream, readyToExportFiles);
		MongoCollection collWithPojoCodec = mongoTemplate.getDb().withCodecRegistry(ExportManager.pojoCodecRegistry).getCollection(tmpVarCollName != null ? tmpVarCollName : mongoTemplate.getCollectionName(VariantRunData.class));

		Map<String, Integer> individualPositions = new LinkedHashMap<>();
		for (String ind : samplesToExport.stream().map(gs -> gs.getIndividual()).distinct().sorted(new AlphaNumericComparator<String>()).collect(Collectors.toList()))
			individualPositions.put(ind, individualPositions.size());

		File warningFile = File.createTempFile("export_warnings_", "");
        FileWriter warningFileWriter = new FileWriter(warningFile);

        String exportName = sModule + "__" + markerCount + "variants__" + individualPositions.size() + "individuals";
        
        if (individualMetadataFieldsToExport != null && !individualMetadataFieldsToExport.isEmpty()) {
        	zos.putNextEntry(new ZipEntry(exportName + ".metadata.tsv"));
        	zos.write("individual".getBytes());
	        IExportHandler.writeMetadataFile(sModule, individualPositions.keySet(), individualMetadataFieldsToExport, zos);
	    	zos.closeEntry();
        }

        zos.putNextEntry(new ZipEntry(exportName + ".gff3"));
        String header = "##gff-version 3" + LINE_SEPARATOR;
        zos.write(header.getBytes());

        TreeMap<String, String> typeToOntology = new TreeMap<>();
        typeToOntology.put(Type.SNP.toString(), "SO:0000694");
        typeToOntology.put(Type.INDEL.toString(), "SO:1000032");
        typeToOntology.put(Type.MIXED.toString(), "SO:0001059");
        typeToOntology.put(Type.SYMBOLIC.toString(), "SO:0000109");
        typeToOntology.put(Type.MNP.toString(), "SO:0001059");

		final Map<Integer, String> sampleIdToIndividualMap = new HashMap<>();
		for (GenotypingSample gs : samplesToExport)
			sampleIdToIndividualMap.put(gs.getId(), gs.getIndividual());
		
        ArrayList<Comparable> unassignedMarkers = new ArrayList<>();
		int nQueryChunkSize = IExportHandler.computeQueryChunkSize(mongoTemplate, markerCount);

//		AtomicLong timeConverting = new AtomicLong(0), timeWriting = new AtomicLong(0);
		
		AbstractExportWritingThread writingThread = new AbstractExportWritingThread() {
			public void run() {
				for (String idOfVarToWrite : markerRunsToWrite.keySet()) {
					if (progress.isAborted() || progress.getError() != null)
						return;

					List<VariantRunData> runsToWrite = markerRunsToWrite.get(idOfVarToWrite);
					if (runsToWrite.isEmpty())
						continue;
					
					List<String> variantDataOrigin = new ArrayList<>();

					StringBuffer sb = new StringBuffer();
					try
					{
						long b4 = System.currentTimeMillis();
						
		                if (markerSynonyms != null) {
		                	String syn = markerSynonyms.get(idOfVarToWrite);
		                    if (syn != null)
		                    	idOfVarToWrite = syn;
		                }

		                VariantRunData vrd = runsToWrite.get(0);	
		                LinkedHashSet<String>[] individualGenotypes = new LinkedHashSet[individualPositions.size()];
	                	for (VariantRunData run : runsToWrite) {
	                    	for (Integer sampleId : run.getSampleGenotypes().keySet()) {
								SampleGenotype sampleGenotype = run.getSampleGenotypes().get(sampleId);
	                            String gtCode = sampleGenotype.getCode();
	                            String individualId = sampleIdToIndividualMap.get(sampleId);
	                            
								if (gtCode == null || !VariantData.gtPassesVcfAnnotationFilters(individualId, sampleGenotype, individuals1, annotationFieldThresholds, individuals2, annotationFieldThresholds2))
									continue;	// skip genotype
								
								int individualIndex = individualPositions.get(individualId);
								if (individualGenotypes[individualIndex] == null)
									individualGenotypes[individualIndex] = new LinkedHashSet<String>();
								individualGenotypes[individualIndex].add(gtCode);
	                        }
	                    }

	                	String refAllele;
						try {
							refAllele = vrd.getKnownAlleleList().get(0);
						}
						catch (IndexOutOfBoundsException ioobe) {	// VariantRunData's known-allele-list is not up to date
							vrd.setKnownAlleleList(mongoTemplate.findById(vrd.getId().getVariantId(), VariantData.class).getKnownAlleleList());
							mongoTemplate.save(vrd);
							refAllele = vrd.getKnownAlleleList().get(0);
						}
		                String chrom = vrd.getReferencePosition() == null ? "0" : vrd.getReferencePosition().getSequence();
		                long start = vrd.getReferencePosition() == null ? 0 : vrd.getReferencePosition().getStartSite();
		                long end = Type.SNP.equals(vrd.getType()) ? start : start + refAllele.length() - 1;
		                sb.append(chrom + "\t" + StringUtils.join(variantDataOrigin, ";") /*source*/ + "\t" + typeToOntology.get(vrd.getType()) + "\t" + start + "\t" + end + "\t" + "." + "\t" + "+" + "\t" + "." + "\t");
		                Comparable syn = markerSynonyms == null ? null : markerSynonyms.get(vrd.getId());
		                sb.append("ID=" + idOfVarToWrite + ";" + (syn != null ? "Name=" + syn + ";" : "") + "alleles=" + StringUtils.join(vrd.getKnownAlleleList(), "/") + ";" + "refallele=" + refAllele + ";");

	            		HashMap<String, List<String>> genotypeStringCache = new HashMap<>();
		                for (String individual : individualPositions.keySet() /* we use this list because it has the proper ordering */) {
		                    int individualIndex = individualPositions.get(individual);
		                    NumberFormat nf = NumberFormat.getInstance(Locale.US);
		                    nf.setMaximumFractionDigits(4);
		                    HashMap<String, Integer> compt = new HashMap<>();
		                    int highestGenotypeCount = 0;
		                    int sum = 0;

		                    HashMap<Object, Integer> genotypeCounts = new HashMap<>(); // will help us to keep track of missing genotypes

		                    String mostFrequentGenotype = null;
		                    if (individualGenotypes[individualIndex] != null) {
		                        for (String genotype : individualGenotypes[individualIndex]) {
		                            if (genotype == null)
		                                continue; /* skip missing genotypes */

		                            int count = 0;
		                            for (String allele : getAllelesFromGenotypeCodeUsingCache(genotype, genotypeStringCache, vrd, mongoTemplate)) {
		                                for (String ka : vrd.getKnownAlleleList()) {
		                                    if (allele.equals(ka) && !(compt.containsKey(ka))) {
		                                        count++;
		                                        compt.put(ka, count);
		                                    } else if (allele.equals(ka) && compt.containsKey(ka)) {
		                                        if (compt.get(ka) != 0) {
		                                            count++;
		                                            compt.put(ka, count);
		                                        } else {
		                                            compt.put(ka, count);
		                                        }
		                                    } else if (!(compt.containsKey(ka))) {
		                                        compt.put(ka, 0);
		                                    }
		                                }
		                            }
		                            for (int countValue : compt.values()) {
		                                sum += countValue;
		                            }

		                            int gtCount = 1 + Helper.getCountForKey(genotypeCounts, genotype);
		                            if (gtCount > highestGenotypeCount) {
		                                highestGenotypeCount = gtCount;
		                                mostFrequentGenotype = genotype;
		                            }
		                            genotypeCounts.put(genotype, gtCount);
		                        }
		                    }

		        			List<String> alleles = getAllelesFromGenotypeCodeUsingCache(mostFrequentGenotype, genotypeStringCache, vrd, mongoTemplate);
		                    if (!alleles.isEmpty()) {
		                        sb.append("acounts=" + individual + ":");

		                        for (String knowAllelesCompt : compt.keySet()) {
		                            sb.append(knowAllelesCompt + " " + nf.format(compt.get(knowAllelesCompt) / (float) sum) + " " + compt.get(knowAllelesCompt) + " ");
		                        }
		                        sb.append(alleles.size() + ";");
		                    }
		                    if (genotypeCounts.size() > 1) {
		                        String sVariantId = markerSynonyms != null ? markerSynonyms.get(idOfVarToWrite) : idOfVarToWrite;
		                        warningFileWriter.write("- Dissimilar genotypes found for variant " + (sVariantId == null ? vrd.getId() : sVariantId) + ", individual " + individual + ". Exporting most frequent: " + StringUtils.join(alleles, ",") + "\n");
		                    }
		                }
		                sb.append(LINE_SEPARATOR);
		                zos.write(sb.toString().getBytes());
		            }
					catch (Exception e)
					{
						if (progress.getError() == null)	// only log this once
							LOG.debug("Unable to export " + idOfVarToWrite, e);
						progress.setError("Unable to export " + idOfVarToWrite + ": " + e.getMessage());
					}
				}
				markerRunsToWrite.clear();
			}
		};
		
		ExportManager exportManager = new ExportManager(mongoTemplate, collWithPojoCodec, VariantRunData.class, varQuery, samplesToExport, true, nQueryChunkSize, writingThread, markerCount, warningFileWriter, progress);
		exportManager.readAndWrite();	
        zos.closeEntry();
        
        if (unassignedMarkers.size() > 0)
        	LOG.info("No chromosomal position found for " + unassignedMarkers.size() + " markers " + StringUtils.join(unassignedMarkers, ", "));

        warningFileWriter.close();
        if (warningFile.length() > 0) {
            zos.putNextEntry(new ZipEntry(exportName + "-REMARKS.txt"));
            int nWarningCount = 0;
            BufferedReader in = new BufferedReader(new FileReader(warningFile));
            String sLine;
            while ((sLine = in.readLine()) != null) {
                zos.write((sLine + "\n").getBytes());
                nWarningCount++;
            }
            LOG.info("Number of Warnings for export (" + exportName + "): " + nWarningCount);
            in.close();
            zos.closeEntry();
        }
        warningFile.delete();

        zos.finish();
        zos.close();
        progress.setCurrentStepProgress((short) 100);
    }

	private List<String> getAllelesFromGenotypeCodeUsingCache(String gtCode, HashMap<String, List<String>> genotypeStringCache, VariantRunData vrd, MongoTemplate mongoTemplate) {
		List<String> alleles = genotypeStringCache.get(gtCode);
        if (alleles == null) {
        	alleles = vrd.safelyGetAllelesFromGenotypeCode(gtCode, mongoTemplate);
        	genotypeStringCache.put(gtCode, alleles);
        }
        return alleles;
	}
	
    /* (non-Javadoc)
 * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
     */
    @Override
    public List<String> getStepList() {
        return Arrays.asList(new String[]{"Exporting data to GFF3 format"});
    }
    
	@Override
	public String[] getExportDataFileExtensions() {
		return new String[] {"gff3"};
	}
}
