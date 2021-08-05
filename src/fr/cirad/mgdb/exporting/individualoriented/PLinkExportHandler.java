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
package fr.cirad.mgdb.exporting.individualoriented;

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.exporting.tools.ExportManager;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongodao.MgdbDao;
import fr.cirad.tools.Helper;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mongo.MongoTemplateManager;
import htsjdk.variant.variantcontext.VariantContext.Type;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.bson.Document;
import org.springframework.data.mongodb.core.MongoTemplate;

import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoCursor;

/**
 * The Class PLinkExportHandler.
 */
public class PLinkExportHandler extends AbstractIndividualOrientedExportHandler {

    /**
     * The Constant LOG.
     */
    private static final Logger LOG = Logger.getLogger(PLinkExportHandler.class);

    /**
     * The supported variant types.
     */
    private static List<String> supportedVariantTypes;

    static {
        supportedVariantTypes = new ArrayList<String>();
        supportedVariantTypes.add(Type.SNP.toString());
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatName()
     */
    @Override
    public String getExportFormatName() {
        return "PLINK";
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
     */
    @Override
    public String getExportFormatDescription() {
    	return "Exports zipped PED and MAP files. See <a target='_blank' href='http://zzz.bwh.harvard.edu/plink/data.shtml#ped'>http://zzz.bwh.harvard.edu/plink/data.shtml#ped</a> for more details";
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.individualoriented.AbstractIndividualOrientedExportHandler#getSupportedVariantTypes()
     */
    @Override
    public List<String> getSupportedVariantTypes() {
        return supportedVariantTypes;
    }
    
	@Override
	public String getExportArchiveExtension() {
		return "zip";
	}

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.individualoriented.AbstractIndividualOrientedExportHandler#exportData(java.io.OutputStream, java.lang.String, java.util.Collection, boolean, fr.cirad.tools.ProgressIndicator, com.mongodb.DBCursor, java.util.Map, java.util.Map)
     */
    @Override
    public void exportData(OutputStream outputStream, String sModule, File[] individualExportFiles, boolean fDeleteSampleExportFilesOnExit, ProgressIndicator progress, String tmpVarCollName, Document varQuery, long markerCount, Map<String, String> markerSynonyms, Collection<String> individualMetadataFieldsToExport, Map<String, InputStream> readyToExportFiles) throws Exception {
		MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        int nQueryChunkSize = IExportHandler.computeQueryChunkSize(mongoTemplate, markerCount);
        File warningFile = File.createTempFile("export_warnings_", "");
        FileWriter warningFileWriter = new FileWriter(warningFile);
        ZipOutputStream zos = IExportHandler.createArchiveOutputStream(outputStream, readyToExportFiles);
		MongoCollection collWithPojoCodec = mongoTemplate.getDb().withCodecRegistry(ExportManager.pojoCodecRegistry).getCollection(tmpVarCollName != null ? tmpVarCollName : mongoTemplate.getCollectionName(VariantRunData.class));
        String exportName = sModule + "__" + markerCount + "variants__" + individualExportFiles.length + "individuals";
        
        if (individualMetadataFieldsToExport != null && !individualMetadataFieldsToExport.isEmpty()) {
        	zos.putNextEntry(new ZipEntry(sModule + "__" + individualExportFiles.length + "individuals_metadata.tsv"));
        	zos.write("individual".getBytes());
	        ArrayList<String> exportedIndividuals = new ArrayList<>();
	        for (File indFile : individualExportFiles)
	        	try (Scanner scanner = new Scanner(indFile)) {
	        		exportedIndividuals.add(scanner.nextLine());
	        	}
	        IExportHandler.writeMetadataFile(sModule, exportedIndividuals, individualMetadataFieldsToExport, zos);
	    	zos.closeEntry();
        }
        
        zos.putNextEntry(new ZipEntry(exportName + ".ped"));
        TreeMap<Integer, Comparable> problematicMarkerIndexToNameMap = writeGenotypeFile(zos, sModule, nQueryChunkSize, null, varQuery, markerSynonyms, individualExportFiles, warningFileWriter, progress);
    	zos.closeEntry();

        zos.putNextEntry(new ZipEntry(exportName + ".map"));

        int nMarkerIndex = 0;
        ArrayList<Comparable> unassignedMarkers = new ArrayList<>();
		try (MongoCursor<Document> markerCursor = IExportHandler.getMarkerCursorWithCorrectCollation(mongoTemplate.getCollection(tmpVarCollName != null ? tmpVarCollName : mongoTemplate.getCollectionName(VariantData.class)), varQuery, nQueryChunkSize)) {
	        while (markerCursor.hasNext()) {
	            Document exportVariant = markerCursor.next();
	            Document refPos = (Document) exportVariant.get(VariantData.FIELDNAME_REFERENCE_POSITION);
	            Long pos = refPos == null ? null : ((Number) refPos.get(ReferencePosition.FIELDNAME_START_SITE)).longValue();
	            String markerId = (String) exportVariant.get("_id");
	            String chrom = refPos == null ? null : (String) refPos.get(ReferencePosition.FIELDNAME_SEQUENCE);
	            if (chrom == null)
	            	unassignedMarkers.add(markerId);
	            String exportedId = markerSynonyms == null ? markerId : markerSynonyms.get(markerId);
	            zos.write(((chrom == null ? "0" : chrom) + " " + exportedId + " " + 0 + " " + (pos == null ? 0 : pos) + LINE_SEPARATOR).getBytes());
	
	            if (problematicMarkerIndexToNameMap.containsKey(nMarkerIndex)) {	// we are going to need this marker's name for the warning file
	            	String variantName = markerId;
	                if (markerSynonyms != null) {
	                	String syn = markerSynonyms.get(markerId);
	                    if (syn != null)
	                        variantName = syn;
	                }
	                problematicMarkerIndexToNameMap.put(nMarkerIndex, variantName);
	            }
	            nMarkerIndex++;
	        }
		}
        zos.closeEntry();
        
        if (unassignedMarkers.size() > 0)
        	LOG.info("No chromosomal position found for " + unassignedMarkers.size() + " markers " + StringUtils.join(unassignedMarkers, ", "));

        if (warningFile.length() > 0) {
            zos.putNextEntry(new ZipEntry(exportName + "-REMARKS.txt"));
            int nWarningCount = 0;
            BufferedReader in = new BufferedReader(new FileReader(warningFile));
            String sLine;
            while ((sLine = in.readLine()) != null) {
                for (Integer aMarkerIndex : problematicMarkerIndexToNameMap.keySet()) {
                    sLine = sLine.replaceAll("__" + aMarkerIndex + "__", problematicMarkerIndexToNameMap.get(aMarkerIndex).toString());
                }
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
    
    public TreeMap<Integer, Comparable> writeGenotypeFile(OutputStream os, String sModule, int nQueryChunkSize, MongoCollection<Document> varColl, Document varQuery, Map<String, String> markerSynonyms, File[] individualExportFiles, FileWriter warningFileWriter, ProgressIndicator progress) throws IOException {
        TreeMap<Integer, Comparable> problematicMarkerIndexToNameMap = new TreeMap<Integer, Comparable>();
        short nProgress = 0, nPreviousProgress = 0;
        int i = 0;
        try
        {
	        for (File f : individualExportFiles) {
	            BufferedReader in = new BufferedReader(new FileReader(f));
	            try {
	                String individualId, line = in.readLine();	// read sample id
	                if (line != null) {
	                    individualId = line;
	                    String population = MgdbDao.getIndividualPopulation(sModule, line);
	                    String individualInfo = (population == null ? "." : population) + " " + individualId;
	                    os.write((individualInfo + " 0 0 0 " + getIndividualGenderCode(sModule, individualId)).getBytes());
	                } else {
	                    throw new Exception("Unable to read first line of temp export file " + f.getName());
	                }
	
	                int nMarkerIndex = 0;
	                while ((line = in.readLine()) != null) {
	                    List<String> genotypes = Helper.split(line, "|");
	                    HashMap<Object, Integer> genotypeCounts = new HashMap<>();	// will help us to keep track of missing genotypes
	                    int highestGenotypeCount = 0;
	                    String mostFrequentGenotype = null;
	                    for (String genotype : genotypes) {
	                        if (genotype.length() == 0)
	                            continue;	/* skip missing genotypes */
	
	                        int gtCount = 1 + Helper.getCountForKey(genotypeCounts, genotype);
	                        if (gtCount > highestGenotypeCount) {
	                            highestGenotypeCount = gtCount;
	                            mostFrequentGenotype = genotype;
	                        }
	                        genotypeCounts.put(genotype, gtCount);
	                    }
	
	                    if (genotypeCounts.size() > 1) {
	                    	if (warningFileWriter != null)
	                    		warningFileWriter.write("- Dissimilar genotypes found for variant " + nMarkerIndex + ", individual " + individualId + ". Exporting most frequent: " + mostFrequentGenotype + "\n");
	                        problematicMarkerIndexToNameMap.put(nMarkerIndex, "");
	                    }
	
	                    String[] alleles = mostFrequentGenotype == null ? new String[0] : mostFrequentGenotype.split(" ");
	                    if (alleles.length > 2) {
	                    	if (warningFileWriter != null)
	                    		warningFileWriter.write("- More than 2 alleles found for variant " + nMarkerIndex + ", individual " + individualId + ". Exporting only the first 2 alleles.\n");
	                        problematicMarkerIndexToNameMap.put(nMarkerIndex, "");
	                    }
	
	                    String all1 = alleles.length == 0 ? "0" : alleles[0];
	                    String all2 = alleles.length == 0 ? "0" : alleles[alleles.length == 1 ? 0 : 1];
	                    if (all1.length() != 1 || all2.length() != 1) {
	                    	if (warningFileWriter != null)
	                    		warningFileWriter.write("- SNP expected, but alleles are not coded on a single char for variant " + nMarkerIndex + ", individual " + individualId + ". Ignoring this genotype.\n");
	                        problematicMarkerIndexToNameMap.put(nMarkerIndex, "");
	                    } else {
	                        os.write((" " + all1 + " " + all2).getBytes());
	                    }
	
	                    nMarkerIndex++;
	                }
	            } catch (Exception e) {
	                LOG.error("Error exporting data", e);
	                progress.setError("Error exporting data: " + e.getClass().getSimpleName() + (e.getMessage() != null ? " - " + e.getMessage() : ""));
	                return null;
	            } finally {
	                in.close();               
	            }
	
	            if (progress.isAborted()) {
	                return null;
	            }
	
	            nProgress = (short) (++i * 100 / individualExportFiles.length);
	            if (nProgress > nPreviousProgress) {
	                progress.setCurrentStepProgress(nProgress);
	                nPreviousProgress = nProgress;
	            }
	            os.write('\n');
	        }
        }
        finally
        {
        	for (File f : individualExportFiles)
                if (!f.delete()) {
                    f.deleteOnExit();
                    LOG.info("Unable to delete tmp export file " + f.getAbsolutePath());
                }
        }
    	if (warningFileWriter != null)
    		warningFileWriter.close();
    	
        return problematicMarkerIndexToNameMap;
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
     */
    @Override
    public List<String> getStepList() {
        return Arrays.asList(new String[]{"Exporting data to PLINK format"});
    }

    /**
     * Gets the individual gender code.
     *
     * @param sModule the module
     * @param individual the individual
     * @return the individual gender code
     */
    protected String getIndividualGenderCode(String sModule, String individual) {
        return "U";
    }
    
	@Override
	public String[] getExportDataFileExtensions() {
		return new String[] {"ped", "map"};
	}
}
