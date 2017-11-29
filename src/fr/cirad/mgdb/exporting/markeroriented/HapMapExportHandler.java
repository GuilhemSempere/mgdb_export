/**
 * *****************************************************************************
 * MGDB Export - Mongo Genotype DataBase, export handlers
 * Copyright (C) 2016 <CIRAD>
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
 * See <http://www.gnu.org/licenses/agpl.html> for details about GNU Affero
 * General Public License V3.
 ******************************************************************************
 */
package fr.cirad.mgdb.exporting.markeroriented;

import fr.cirad.mgdb.model.mongo.maintypes.Individual;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongo.subtypes.SampleGenotype;
import fr.cirad.mgdb.model.mongo.subtypes.SampleId;
import fr.cirad.mgdb.model.mongodao.MgdbDao;
import fr.cirad.tools.AlphaNumericComparator;
import fr.cirad.tools.Helper;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mongo.MongoTemplateManager;
import htsjdk.variant.variantcontext.VariantContext.Type;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.springframework.data.mongodb.core.MongoTemplate;

import com.mongodb.DBCursor;
import com.mongodb.DBObject;

// TODO: Auto-generated Javadoc
/**
 * The Class HapMapExportHandler.
 */
public class HapMapExportHandler extends AbstractMarkerOrientedExportHandler {

    /**
     * The Constant LOG.
     */
    private static final Logger LOG = Logger.getLogger(HapMapExportHandler.class);

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
        return "HAPMAP";
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
     */
    @Override
    public String getExportFormatDescription() {
        return "Exports data in HapMap Format. See <a target='_blank' href='http://heidi.chnebu.ch/doku.php?id=hapmap'>http://heidi.chnebu.ch/doku.php?id=hapmap</a> for more details";
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.markeroriented.AbstractMarkerOrientedExportHandler#getSupportedVariantTypes()
     */
    @Override
    public List<String> getSupportedVariantTypes() {
        return supportedVariantTypes;
    }
    
	@Override
	public String getExportFileExtension() {
		return "zip";
	}

    /* (non-Javadoc)
 * @see fr.cirad.mgdb.exporting.markeroriented.AbstractMarkerOrientedExportHandler#exportData(java.io.OutputStream, java.lang.String, java.util.List, fr.cirad.tools.ProgressIndicator, com.mongodb.DBCursor, java.util.Map, int, int, java.util.Map)
     */
    @Override
    public void exportData(OutputStream outputStream, String sModule, List<SampleId> sampleIDs1, List<SampleId> sampleIDs2, ProgressIndicator progress, DBCursor markerCursor, Map<Comparable, Comparable> markerSynonyms, HashMap<String, Integer> annotationFieldThresholds, HashMap<String, Integer> annotationFieldThresholds2, Map<String, InputStream> readyToExportFiles) throws Exception {
        
		List<String> individuals1 = getIndividualsFromSamples(sModule, sampleIDs1).stream().map(ind -> ind.getId()).collect(Collectors.toList());	
		List<String> individuals2 = getIndividualsFromSamples(sModule, sampleIDs2).stream().map(ind -> ind.getId()).collect(Collectors.toList());

		ArrayList<SampleId> sampleIDs = (ArrayList<SampleId>) CollectionUtils.union(sampleIDs1, sampleIDs2);
		List<Individual> individuals = getIndividualsFromSamples(sModule, sampleIDs), sortedIndividuals = new ArrayList<Individual>(individuals);
		Collections.sort(sortedIndividuals, new AlphaNumericComparator<Individual>());

		MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        File warningFile = File.createTempFile("export_warnings_", "");
        FileWriter warningFileWriter = new FileWriter(warningFile);

        int markerCount = markerCursor.count();

        ZipOutputStream zos = new ZipOutputStream(outputStream);

        if (readyToExportFiles != null) {
            for (String readyToExportFile : readyToExportFiles.keySet()) {
                zos.putNextEntry(new ZipEntry(readyToExportFile));
                InputStream inputStream = readyToExportFiles.get(readyToExportFile);
                byte[] dataBlock = new byte[1024];
                int count = inputStream.read(dataBlock, 0, 1024);
                while (count != -1) {
                    zos.write(dataBlock, 0, count);
                    count = inputStream.read(dataBlock, 0, 1024);
                }
                zos.closeEntry();
            }
        }

        String exportName = sModule + "_" + markerCount + "variants_" + individuals.size() + "individuals";
        zos.putNextEntry(new ZipEntry(exportName + ".hapmap"));
        String header = "rs#" + "\t" + "alleles" + "\t" + "chrom" + "\t" + "pos" + "\t" + "strand" + "\t" + "assembly#" + "\t" + "center" + "\t" + "protLSID" + "\t" + "assayLSID" + "\t" + "panelLSID" + "\t" + "QCcode";
        zos.write(header.getBytes());
        for (int i = 0; i < sortedIndividuals.size(); i++) {
            zos.write(("\t" + sortedIndividuals.get(i).getId()).getBytes());
        }
        zos.write((LINE_SEPARATOR).getBytes());

        int avgObjSize = (Integer) mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantRunData.class)).getStats().get("avgObjSize");
        int nChunkSize = nMaxChunkSizeInMb * 1024 * 1024 / avgObjSize;
        short nProgress = 0, nPreviousProgress = 0;
        long nLoadedMarkerCount = 0;

        while (markerCursor == null || markerCursor.hasNext()) {
            int nLoadedMarkerCountInLoop = 0;
            Map<Comparable, String> markerChromosomalPositions = new LinkedHashMap<Comparable, String>();
            boolean fStartingNewChunk = true;
            markerCursor.batchSize(nChunkSize);
            while (markerCursor.hasNext() && (fStartingNewChunk || nLoadedMarkerCountInLoop % nChunkSize != 0)) {
                DBObject exportVariant = markerCursor.next();
                DBObject refPos = (DBObject) exportVariant.get(VariantData.FIELDNAME_REFERENCE_POSITION);
                markerChromosomalPositions.put((Comparable) exportVariant.get("_id"), refPos == null ? null : (refPos.get(ReferencePosition.FIELDNAME_SEQUENCE) + ":" + refPos.get(ReferencePosition.FIELDNAME_START_SITE)));
                nLoadedMarkerCountInLoop++;
                fStartingNewChunk = false;
            }

            List<Comparable> currentMarkers = new ArrayList<Comparable>(markerChromosomalPositions.keySet());
            LinkedHashMap<VariantData, Collection<VariantRunData>> variantsAndRuns = MgdbDao.getSampleGenotypes(mongoTemplate, sampleIDs, currentMarkers, true, null /*new Sort(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ChromosomalPosition.FIELDNAME_SEQUENCE).and(new Sort(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ChromosomalPosition.FIELDNAME_START_SITE))*/);	// query mongo db for matching genotypes			
            for (VariantData variant : variantsAndRuns.keySet()) // read data and write results into temporary files (one per sample)
            {
                Comparable variantId = variant.getId();
                if (markerSynonyms != null) {
                    Comparable syn = markerSynonyms.get(variantId);
                    if (syn != null) {
                        variantId = syn;
                    }
                }

                boolean fIsSNP = variant.getType().equals(Type.SNP.toString());
                byte[] missingGenotype = ("\t" + "NN").getBytes();

                String refPos = markerChromosomalPositions.get(variant.getId());
                String[] chromAndPos = refPos == null ? null : refPos.split(":");
                zos.write(((variantId == null ? variant.getId() : variantId) + "\t" + StringUtils.join(variant.getKnownAlleleList(), "/") + "\t" + (chromAndPos == null ? 0 : chromAndPos[0]) + "\t" + (chromAndPos == null ? 0 : Long.parseLong(chromAndPos[1])) + "\t" + "+").getBytes());
                for (int j = 0; j < 6; j++) {
                    zos.write(("\t" + "NA").getBytes());
                }

                Map<String, List<String>> individualGenotypes = new TreeMap<String, List<String>>(new AlphaNumericComparator<String>());
                Collection<VariantRunData> runs = variantsAndRuns.get(variant);
                if (runs != null) {
                    for (VariantRunData run : runs) {
                        for (Integer sampleIndex : run.getSampleGenotypes().keySet()) {
                            SampleGenotype sampleGenotype = run.getSampleGenotypes().get(sampleIndex);
                            String individualId = individuals.get(sampleIDs.indexOf(new SampleId(run.getId().getProjectId(), sampleIndex))).getId();
                            
							if (!VariantData.gtPassesVcfAnnotationFilters(individualId, sampleGenotype, individuals1, annotationFieldThresholds, individuals2, annotationFieldThresholds2))
								continue;	// skip genotype
							
                            String gtCode = run.getSampleGenotypes().get(sampleIndex).getCode();
                            List<String> storedIndividualGenotypes = individualGenotypes.get(individualId);
                            if (storedIndividualGenotypes == null) {
                                storedIndividualGenotypes = new ArrayList<String>();
                                individualGenotypes.put(individualId, storedIndividualGenotypes);
                            }
                            storedIndividualGenotypes.add(gtCode);
                        }
                    }
                }

                int writtenGenotypeCount = 0;
                for (Individual individual : sortedIndividuals /* we use this list because it has the proper ordering */) {
                    int individualIndex = individuals.indexOf(individual.getId());
                    while (writtenGenotypeCount < individualIndex - 1) {
                        zos.write(missingGenotype);
                        writtenGenotypeCount++;
                    }

                    List<String> genotypes = individualGenotypes.get(individual.getId());
                    HashMap<Object, Integer> genotypeCounts = new HashMap<Object, Integer>();	// will help us to keep track of missing genotypes
                    int highestGenotypeCount = 0;
                    String mostFrequentGenotype = null;
                    if (genotypes != null) {
                        for (String genotype : genotypes) {
                            if (genotype.length() == 0) {
                                continue;	/* skip missing genotypes */
                            }

                            int gtCount = 1 + Helper.getCountForKey(genotypeCounts, genotype);
                            if (gtCount > highestGenotypeCount) {
                                highestGenotypeCount = gtCount;
                                mostFrequentGenotype = genotype;
                            }
                            genotypeCounts.put(genotype, gtCount);
                        }
                    }

                    byte[] exportedGT = mostFrequentGenotype == null ? missingGenotype : ("\t" + StringUtils.join(variant.getAllelesFromGenotypeCode(mostFrequentGenotype), fIsSNP ? "" : "/")).getBytes();
                    zos.write(exportedGT);
                    writtenGenotypeCount++;

                    if (genotypeCounts.size() > 1) {
                        warningFileWriter.write("- Dissimilar genotypes found for variant " + (variantId == null ? variant.getId() : variantId) + ", individual " + individual.getId() + ". Exporting most frequent: " + new String(exportedGT) + "\n");
                    }
                }

                while (writtenGenotypeCount < individuals.size()) {
                    zos.write(missingGenotype);
                    writtenGenotypeCount++;
                }
                zos.write((LINE_SEPARATOR).getBytes());
            }

            if (progress.hasAborted()) {
                warningFileWriter.close();
                return;
            }

            nLoadedMarkerCount += nLoadedMarkerCountInLoop;
            nProgress = (short) (nLoadedMarkerCount * 100 / markerCount);
            if (nProgress > nPreviousProgress) {
                //				if (nProgress%5 == 0)
                //					LOG.info("========================= exportData: " + nProgress + "% =========================" + (System.currentTimeMillis() - before)/1000 + "s");
                progress.setCurrentStepProgress(nProgress);
                nPreviousProgress = nProgress;
            }
        }
        zos.closeEntry();

        warningFileWriter.close();
        if (warningFile.length() > 0) {
            zos.putNextEntry(new ZipEntry(exportName + "-REMARKS.txt"));
            int nWarningCount = 0;
            BufferedReader in = new BufferedReader(new FileReader(warningFile));
            String sLine;
            while ((sLine = in.readLine()) != null) {
                zos.write((sLine + "\n").getBytes());
                in.readLine();
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

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
     */
    @Override
    public List<String> getStepList() {
        return Arrays.asList(new String[]{"Exporting data to HAPMAP format"});
    }
}
