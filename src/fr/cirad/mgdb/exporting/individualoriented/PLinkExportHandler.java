/**
 * *****************************************************************************
 * MGDB Export - Mongo Genotype DataBase, export handlers
 * Copyright (C) 2016 <CIRAD>
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License, version 3 as published by the
 * Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
 * details.
 *
 * See <http://www.gnu.org/licenses/agpl.html> for details about GNU General
 * Public License V3.
 ******************************************************************************
 */
package fr.cirad.mgdb.exporting.individualoriented;

import fr.cirad.mgdb.model.mongo.maintypes.Individual;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongodao.MgdbDao;
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
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.log4j.Logger;
import org.springframework.data.mongodb.core.MongoTemplate;

import com.mongodb.DBCursor;
import com.mongodb.DBObject;

// TODO: Auto-generated Javadoc
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
        return "Exports zipped PED and MAP files. See <a target='_blank' href='http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml'>http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml</a> for more details";
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.individualoriented.AbstractIndividualOrientedExportHandler#getSupportedVariantTypes()
     */
    @Override
    public List<String> getSupportedVariantTypes() {
        return supportedVariantTypes;
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.individualoriented.AbstractIndividualOrientedExportHandler#exportData(java.io.OutputStream, java.lang.String, java.util.Collection, boolean, fr.cirad.tools.ProgressIndicator, com.mongodb.DBCursor, java.util.Map, java.util.Map)
     */
    @Override
    public void exportData(OutputStream outputStream, String sModule, Collection<File> individualExportFiles, boolean fDeleteSampleExportFilesOnExit, ProgressIndicator progress, DBCursor markerCursor, Map<Comparable, Comparable> markerSynonyms, Map<String, InputStream> readyToExportFiles) throws Exception {
        File warningFile = File.createTempFile("export_warnings_", "");
        FileWriter warningFileWriter = new FileWriter(warningFile);

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

        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        int markerCount = markerCursor.count();

        String exportName = sModule + "_" + markerCount + "variants_" + individualExportFiles.size() + "individuals";
        zos.putNextEntry(new ZipEntry(exportName + ".ped"));

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
	                    String population = getIndividualPopulation(sModule, line);
	                    String individualInfo = (population == null ? "." : population) + " " + individualId;
	                    zos.write((individualInfo + " 0 0 0 " + getIndividualGenderCode(sModule, individualId)).getBytes());
	                } else {
	                    throw new Exception("Unable to read first line of temp export file " + f.getName());
	                }
	
	                int nMarkerIndex = 0;
	                while ((line = in.readLine()) != null) {
	                    List<String> genotypes = MgdbDao.split(line, "|");
	                    HashMap<Object, Integer> genotypeCounts = new HashMap<Object, Integer>();	// will help us to keep track of missing genotypes
	                    int highestGenotypeCount = 0;
	                    String mostFrequentGenotype = null;
	                    for (String genotype : genotypes) {
	                        if (genotype.length() == 0) {
	                            continue;	/* skip missing genotypes */
	                        }
	
	                        int gtCount = 1 + MgdbDao.getCountForKey(genotypeCounts, genotype);
	                        if (gtCount > highestGenotypeCount) {
	                            highestGenotypeCount = gtCount;
	                            mostFrequentGenotype = genotype;
	                        }
	                        genotypeCounts.put(genotype, gtCount);
	                    }
	
	                    if (genotypeCounts.size() > 1) {
	                        warningFileWriter.write("- Dissimilar genotypes found for variant " + nMarkerIndex + ", individual " + individualId + ". Exporting most frequent: " + mostFrequentGenotype + "\n");
	                        problematicMarkerIndexToNameMap.put(nMarkerIndex, "");
	                    }
	
	                    String[] alleles = mostFrequentGenotype == null ? new String[0] : mostFrequentGenotype.split(" ");
	                    if (alleles.length > 2) {
	                        warningFileWriter.write("- More than 2 alleles found for variant " + nMarkerIndex + ", individual " + individualId + ". Exporting only the first 2 alleles.\n");
	                        problematicMarkerIndexToNameMap.put(nMarkerIndex, "");
	                    }
	
	                    String all1 = alleles.length == 0 ? "0" : alleles[0];
	                    String all2 = alleles.length == 0 ? "0" : alleles[alleles.length == 1 ? 0 : 1];
	                    if (all1.length() != 1 || all2.length() != 1) {
	                        warningFileWriter.write("- SNP expected, but alleles are not coded on a single char for variant " + nMarkerIndex + ", individual " + individualId + ". Ignoring this genotype.\n");
	                        problematicMarkerIndexToNameMap.put(nMarkerIndex, "");
	                    } else {
	                        zos.write((" " + all1 + " " + all2).getBytes());
	                    }
	
	                    nMarkerIndex++;
	                }
	            } catch (Exception e) {
	                LOG.error("Error exporting data", e);
	                progress.setError("Error exporting data: " + e.getClass().getSimpleName() + (e.getMessage() != null ? " - " + e.getMessage() : ""));
	                return;
	            } finally {
	                in.close();               
	            }
	
	            if (progress.hasAborted()) {
	                return;
	            }
	
	            nProgress = (short) (++i * 100 / individualExportFiles.size());
	            if (nProgress > nPreviousProgress) {
	                progress.setCurrentStepProgress(nProgress);
	                nPreviousProgress = nProgress;
	            }
	            zos.write('\n');
	        }
        }
        finally
        {
        	for (File f : individualExportFiles)
                if (!f.delete()) {
                    f.deleteOnExit();
                    LOG.info("Unable to delete tmp export file " + f.getAbsolutePath());
                }
        	zos.closeEntry();
        }
        warningFileWriter.close();

        zos.putNextEntry(new ZipEntry(exportName + ".map"));

        int avgObjSize = (Integer) mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantRunData.class)).getStats().get("avgObjSize");
        int nChunkSize = nMaxChunkSizeInMb * 1024 * 1024 / avgObjSize;

        markerCursor.batchSize(nChunkSize);
        int nMarkerIndex = 0;
        while (markerCursor.hasNext()) {
            DBObject exportVariant = markerCursor.next();
            DBObject refPos = (DBObject) exportVariant.get(VariantData.FIELDNAME_REFERENCE_POSITION);
            Comparable markerId = (Comparable) exportVariant.get("_id");
            String chrom = (String) refPos.get(ReferencePosition.FIELDNAME_SEQUENCE);
            Long pos = ((Number) refPos.get(ReferencePosition.FIELDNAME_START_SITE)).longValue();

            if (chrom == null) {
                LOG.warn("Chromosomal position not found for marker " + markerId);
            }
            Comparable exportedId = markerSynonyms == null ? markerId : markerSynonyms.get(markerId);
            zos.write(((chrom == null ? "0" : chrom) + " " + exportedId + " " + 0 + " " + (pos == null ? 0 : pos) + LINE_SEPARATOR).getBytes());

            if (problematicMarkerIndexToNameMap.containsKey(nMarkerIndex)) {	// we are going to need this marker's name for the warning file
                Comparable variantName = markerId;
                if (markerSynonyms != null) {
                    Comparable syn = markerSynonyms.get(markerId);
                    if (syn != null) {
                        variantName = syn;
                    }
                }
                problematicMarkerIndexToNameMap.put(nMarkerIndex, variantName);
            }
            nMarkerIndex++;
        }
        zos.closeEntry();

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
                sLine = in.readLine();
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
        return Arrays.asList(new String[]{"Exporting data to PLINK format"});
    }

    /**
     * Gets the individual population.
     *
     * @param sModule the module
     * @param individual the individual
     * @return the individual population
     */
    protected String getIndividualPopulation(final String sModule, final String individual) {
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        return mongoTemplate.findById(individual, Individual.class).getPopulation();
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
}
