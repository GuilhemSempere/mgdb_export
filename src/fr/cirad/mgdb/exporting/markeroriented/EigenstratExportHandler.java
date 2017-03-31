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
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.springframework.data.mongodb.core.MongoTemplate;

import com.mongodb.DBCursor;
import com.mongodb.DBObject;

// TODO: Auto-generated Javadoc
/**
 * The Class EigenstratExportHandler.
 */
public class EigenstratExportHandler extends AbstractMarkerOrientedExportHandler {

    /**
     * The Constant LOG.
     */
    private static final Logger LOG = Logger.getLogger(EigenstratExportHandler.class);

    /**
     * The Constant EIGENSTRAT_FORMAT.
     */
    public static final String EIGENSTRAT_FORMAT = "EIGENSTRAT";

    /**
     * The supported variant types.
     */
    private static List<String> supportedVariantTypes;

    static {
        supportedVariantTypes = new ArrayList<String>();
        supportedVariantTypes.add(Type.SNP.toString());
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

    /**
     * Gets the populations from samples.
     *
     * @param sModule the module
     * @param sampleIDs the sample i ds
     * @return the populations from samples
     */
    @SuppressWarnings("unchecked")
    protected List<String> getPopulationsFromSamples(final String sModule, final List<SampleId> sampleIDs) {
        ArrayList<String> result = new ArrayList<String>();
        for (Individual individual : getIndividualsFromSamples(sModule, sampleIDs)) {
            result.add(individual.getPopulation());
        }
        return result;
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatName()
     */
    @Override
    public String getExportFormatName() {
        return EIGENSTRAT_FORMAT;
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
     */
    @Override
    public String getExportFormatDescription() {
        return "Exports zipped ind, snp and eigenstratgeno files, along with an optional remark-file. See <a target='_blank' href='https://github.com/argriffing/eigensoft/blob/master/CONVERTF/README'>https://github.com/argriffing/eigensoft/blob/master/CONVERTF/README</a> for more details";
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
	
    // public static void logMemUsage(Logger log)
    // {
    // long maxMem = Runtime.getRuntime().maxMemory()/(1024*1024);
    // long freeMem = Runtime.getRuntime().freeMemory()/(1024*1024);
    // long totalMem = Runtime.getRuntime().totalMemory()/(1024*1024);
    // log.debug("total: " + totalMem + ", free: " + freeMem + ", max: " +
    // maxMem);
    // }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.markeroriented.AbstractMarkerOrientedExportHandler#exportData(java.io.OutputStream, java.lang.String, java.util.List, fr.cirad.tools.ProgressIndicator, com.mongodb.DBCursor, java.util.Map, int, int, java.util.Map)
     */
    @Override
    public void exportData(OutputStream outputStream, String sModule, List<SampleId> sampleIDs, ProgressIndicator progress, DBCursor markerCursor, Map<Comparable, Comparable> markerSynonyms, int nMinimumGenotypeQuality, int nMinimumReadDepth, Map<String, InputStream> readyToExportFiles) throws Exception {
        // long before = System.currentTimeMillis();

        File warningFile = File.createTempFile("export_warnings_", "");
        FileWriter warningFileWriter = new FileWriter(warningFile);
        File snpFile = null;
        FileWriter snpFileWriter = null; 

        try {
            snpFile = File.createTempFile("snpFile", "");
            snpFileWriter = new FileWriter(snpFile);

            ZipOutputStream zos = new ZipOutputStream(outputStream);
//            if (ByteArrayOutputStream.class.isAssignableFrom(outputStream.getClass())) {
//                zos.setLevel(ZipOutputStream.STORED);
//            }

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

            MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
            int markerCount = markerCursor.count();

            List<Individual> individuals = getIndividualsFromSamples(sModule, sampleIDs);

            ArrayList<String> individualList = new ArrayList<String>();
            StringBuffer indFileContents = new StringBuffer();

            for (int i = 0; i < sampleIDs.size(); i++) {
                Individual individual = individuals.get(i);
                if (!individualList.contains(individual.getId())) {
                    individualList.add(individual.getId());
                    indFileContents.append(individual.getId() + "\t" + getIndividualGenderCode(sModule, individual.getId()) + "\t" + (individual.getPopulation() == null ? "." : individual.getPopulation()) + LINE_SEPARATOR);
                }
            }

            String exportName = sModule + "_" + markerCount + "variants_" + individualList.size() + "individuals";
            zos.putNextEntry(new ZipEntry(exportName + ".ind"));
            zos.write(indFileContents.toString().getBytes());
            zos.closeEntry();
            
            zos.putNextEntry(new ZipEntry(exportName + ".eigenstratgeno"));

            int avgObjSize = (Integer) mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantRunData.class)).getStats().get("avgObjSize");
            int nChunkSize = nMaxChunkSizeInMb * 1024 * 1024 / avgObjSize;
            short nProgress = 0, nPreviousProgress = 0;
            long nLoadedMarkerCount = 0;

            ArrayList<Comparable> unassignedMarkers = new ArrayList<>();
            while (markerCursor.hasNext()) {
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

                    List<String> chromAndPos = Helper.split(markerChromosomalPositions.get(variantId), ":");
                    if (chromAndPos.size() == 0)
                    	unassignedMarkers.add(variantId);
                    // LOG.debug(marker + "\t" + (chromAndPos.length == 0 ? "0" : chromAndPos[0]) + "\t" + 0 + "\t" + (chromAndPos.length == 0 ? 0l : Long.parseLong(chromAndPos[1])) + LINE_SEPARATOR);
                    if (markerSynonyms != null) {
                        Comparable syn = markerSynonyms.get(variantId);
                        if (syn != null) {
                            variantId = syn;
                        }
                    }
                    snpFileWriter.write(variantId + "\t" + (chromAndPos.size() == 0 ? "0" : chromAndPos.get(0)) + "\t" + 0 + "\t" + (chromAndPos.size() == 0 ? 0l : Long.parseLong(chromAndPos.get(1))) + LINE_SEPARATOR);

                    Map<String, List<String>> individualGenotypes = new LinkedHashMap<String, List<String>>();
                    Collection<VariantRunData> runs = variantsAndRuns.get(variant);
                    if (runs != null) {
                        for (VariantRunData run : runs) {
                            for (Integer sampleIndex : run.getSampleGenotypes().keySet()) {
                                SampleGenotype sampleGenotype = run.getSampleGenotypes().get(sampleIndex);
                                String individualId = individuals.get(sampleIDs.indexOf(new SampleId(run.getId().getProjectId(), sampleIndex))).getId();

                                Integer gq = null;
                                try {
                                    gq = (Integer) sampleGenotype.getAdditionalInfo().get(VariantData.GT_FIELD_GQ);
                                } catch (Exception ignored) {
                                }
                                if (gq != null && gq < nMinimumGenotypeQuality) {
                                    continue;
                                }

                                Integer dp = null;
                                try {
                                    dp = (Integer) sampleGenotype.getAdditionalInfo().get(VariantData.GT_FIELD_DP);
                                } catch (Exception ignored) {
                                }
                                if (dp != null && dp < nMinimumReadDepth) {
                                    continue;
                                }

                                String gtCode = sampleGenotype.getCode();
                                List<String> storedIndividualGenotypes = individualGenotypes.get(individualId);
                                if (storedIndividualGenotypes == null) {
                                    storedIndividualGenotypes = new ArrayList<String>();
                                    individualGenotypes.put(individualId, storedIndividualGenotypes);
                                }
                                storedIndividualGenotypes.add(gtCode);
                            }
                        }
                    }

                    for (int j = 0; j < individualList.size(); j++ /* we use this list because it has the proper ordering*/) {
                        String individualId = individualList.get(j);
                        List<String> genotypes = individualGenotypes.get(individualId);
                        HashMap<Object, Integer> genotypeCounts = new HashMap<Object, Integer>(); // will help us to keep track of missing genotypes
                        int highestGenotypeCount = 0;
                        String mostFrequentGenotype = null;
                        if (genotypes != null) {
                            for (String genotype : genotypes) {
                                if (genotype.length() == 0) {
                                    continue; /* skip missing genotypes */
                                }

                                int gtCount = 1 + Helper.getCountForKey(genotypeCounts, genotype);
                                if (gtCount > highestGenotypeCount) {
                                    highestGenotypeCount = gtCount;
                                    mostFrequentGenotype = genotype;
                                }
                                genotypeCounts.put(genotype, gtCount);
                            }
                        }

                        List<String> alleles = mostFrequentGenotype == null ? new ArrayList<String>() : variant.getAllelesFromGenotypeCode(mostFrequentGenotype);

                        int nOutputCode = 0;
                        if (mostFrequentGenotype == null) {
                            nOutputCode = 9;
                        } else {
                            for (String all : Helper.split(mostFrequentGenotype, "/")) {
                                if ("0".equals(all)) {
                                    nOutputCode++;
                                }
                            }
                        }
                        if (j == 0 && variant.getKnownAlleleList().size() > 2) {
                            warningFileWriter.write("- Variant " + variant.getId() + " is multi-allelic. Make sure Eigenstrat genotype encoding specifications are suitable for you.\n");
                        }
                        zos.write(("" + nOutputCode).getBytes());

                        if (genotypeCounts.size() > 1 || alleles.size() > 2) {
                            if (genotypeCounts.size() > 1) {
                                warningFileWriter.write("- Dissimilar genotypes found for variant " + (variantId == null ? variant.getId() : variantId) + ", individual " + individualId + ". Exporting most frequent: " + nOutputCode + "\n");
                            }
                            if (alleles.size() > 2) {
                                warningFileWriter.write("- More than 2 alleles found for variant " + (variantId == null ? variant.getId() : variantId) + ", individual " + individualId + ". Exporting only the first 2 alleles.\n");
                            }
                        }
                    }
                    zos.write((LINE_SEPARATOR).getBytes());
                }

                if (progress.hasAborted()) {
                    warningFileWriter.close();
                    snpFileWriter.close();
                    return;
                }

                nLoadedMarkerCount += nLoadedMarkerCountInLoop;
                nProgress = (short) (nLoadedMarkerCount * 100 / markerCount);
                if (nProgress > nPreviousProgress) {
                    // if (nProgress%5 == 0)
                    // 	LOG.info("============= exportData: " + nProgress + "% =============" + (System.currentTimeMillis() - before)/1000 + "s");
                    progress.setCurrentStepProgress(nProgress);
                    nPreviousProgress = nProgress;
                }
            }
            zos.closeEntry();
            snpFileWriter.close();
            
            if (unassignedMarkers.size() > 0)
            	LOG.info("No chromosomal position found for " + unassignedMarkers.size() + " markers " + StringUtils.join(unassignedMarkers, ", "));
            
            zos.putNextEntry(new ZipEntry(exportName + ".snp"));
            BufferedReader in = new BufferedReader(new FileReader(snpFile));
            String sLine;
            while ((sLine = in.readLine()) != null) {
                zos.write((sLine + "\n").getBytes());
            }
            in.close();
            zos.closeEntry();

            warningFileWriter.close();
            if (warningFile.length() > 0) {
                zos.putNextEntry(new ZipEntry(exportName + "-REMARKS.txt"));
                int nWarningCount = 0;
                in = new BufferedReader(new FileReader(warningFile));
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
        } finally {
            if (snpFile != null && snpFile.exists()) {
                snpFile.delete();
            }
        }
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
     */
    @Override
    public List<String> getStepList() {
        return Arrays.asList(new String[]{"Exporting data to EIGENSTRAT format"});
    }
}
