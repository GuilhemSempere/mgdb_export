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
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;

import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoCursor;

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.model.mongo.maintypes.Assembly;
import fr.cirad.mgdb.model.mongo.maintypes.Individual;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData.VariantRunDataId;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.tools.Helper;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mongo.MongoTemplateManager;

/**
 * The Class FlapjackExportHandler.
 */
public class FlapjackExportHandler extends AbstractIndividualOrientedExportHandler {

    /**
     * The Constant LOG.
     */
    private static final Logger LOG = Logger.getLogger(FlapjackExportHandler.class);

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatName()
     */
    @Override
    public String getExportFormatName() {
        return "FLAPJACK";
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
     */
    @Override
    public String getExportFormatDescription() {
        return "Exports zipped GENOTYPE and MAP files. See <a target='_blank' href='https://ics.hutton.ac.uk/wiki/index.php/Flapjack_Help_-_Projects_and_Data_Formats'>https://ics.hutton.ac.uk/wiki/index.php/Flapjack_Help_-_Projects_and_Data_Formats</a> for more details";
    }
    
	@Override
	public String getExportArchiveExtension() {
		return "fjzip";
	}

    @Override
    public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId, File[] individualExportFiles, boolean fDeleteSampleExportFilesOnExit, ProgressIndicator progress, String tmpVarCollName, Document varQuery, long markerCount, Map<String, String> markerSynonyms, Collection<String> individualMetadataFieldsToExport, Map<String, InputStream> readyToExportFiles) throws Exception {
        File warningFile = File.createTempFile("export_warnings_", "");
        FileWriter warningFileWriter = new FileWriter(warningFile);

        ZipOutputStream zos = IExportHandler.createArchiveOutputStream(outputStream, readyToExportFiles);
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        
        Assembly assembly = mongoTemplate.findOne(new Query(Criteria.where("_id").is(nAssemblyId)), Assembly.class);
        String exportName = sModule + "__" + assembly.getName() + "__" + markerCount + "variants__" + individualExportFiles.length + "individuals";
        int nQueryChunkSize = IExportHandler.computeQueryChunkSize(mongoTemplate, markerCount);
        MongoCollection<Document> varColl = mongoTemplate.getCollection(tmpVarCollName != null ? tmpVarCollName : mongoTemplate.getCollectionName(VariantData.class));

        if (individualMetadataFieldsToExport != null && !individualMetadataFieldsToExport.isEmpty()) {
	        zos.putNextEntry(new ZipEntry(exportName + ".phenotype"));
	        zos.write(("# fjFile = PHENOTYPE" + LINE_SEPARATOR).getBytes());
	        ArrayList<String> exportedIndividuals = new ArrayList<>();
	        for (File indFile : individualExportFiles)
	        	try (Scanner scanner = new Scanner(indFile)) {
	        		exportedIndividuals.add(scanner.nextLine());
	        	}
	        IExportHandler.writeMetadataFile(sModule, exportedIndividuals, individualMetadataFieldsToExport, zos);
	    	zos.closeEntry();
        }

        zos.putNextEntry(new ZipEntry(exportName + ".genotype"));
        TreeMap<Integer, Comparable> problematicMarkerIndexToNameMap = writeGenotypeFile(zos, sModule, nAssemblyId, nQueryChunkSize, varColl, varQuery, markerSynonyms, individualExportFiles, warningFileWriter, progress);
    	zos.closeEntry();

        zos.putNextEntry(new ZipEntry(exportName + ".map"));
        zos.write(("# fjFile = MAP" + LINE_SEPARATOR).getBytes());
        int nMarkerIndex = 0;
        ArrayList<String> unassignedMarkers = new ArrayList<>();
        try (MongoCursor<Document> markerCursor = IExportHandler.getMarkerCursorWithCorrectCollation(varColl, varQuery, nAssemblyId, nQueryChunkSize)) {
            progress.addStep("Writing map file");
            progress.moveToNextStep();
	        while (markerCursor.hasNext()) {
	            Document exportVariant = markerCursor.next();
	            Document refPos = (Document) Helper.readPossiblyNestedField(exportVariant, VariantData.FIELDNAME_REFERENCE_POSITION + "." + nAssemblyId + "." + ReferencePosition.FIELDNAME_SEQUENCE, ";", null);
	            String chrom = refPos == null ? null : (String) refPos.get(ReferencePosition.FIELDNAME_SEQUENCE);
	            Long pos = refPos == null ? null : ((Number) refPos.get(ReferencePosition.FIELDNAME_START_SITE)).longValue();
                String markerId = (String) exportVariant.get("_id");
	            if (chrom == null) 
	            	unassignedMarkers.add(markerId);
	
	            String exportedId = markerSynonyms == null ? markerId : markerSynonyms.get(markerId);
	            zos.write((exportedId + "\t" + (chrom == null ? "0" : chrom) + "\t" + (pos == null ? 0 : pos) + LINE_SEPARATOR).getBytes());
	
//	            if (problematicMarkerIndexToNameMap.containsKey(nMarkerIndex)) {	// we are going to need this marker's name for the warning file
//	            	String variantName = markerId;
//	                if (markerSynonyms != null) {
//	                	String syn = markerSynonyms.get(markerId);
//	                    if (syn != null)
//	                        variantName = syn;
//	                }
//	                problematicMarkerIndexToNameMap.put(nMarkerIndex, variantName);
//	            }
	            progress.setCurrentStepProgress(nMarkerIndex++ * 100 / markerCount);
	        }
	    }
        zos.closeEntry();
        
        if (unassignedMarkers.size() > 0)
        	LOG.info("No chromosomal position found for " + unassignedMarkers.size() + " markers " + StringUtils.join(unassignedMarkers, ", "));

        if (warningFile.length() > 0) {
            progress.addStep("Adding lines to warning file");
            progress.moveToNextStep();
            progress.setPercentageEnabled(false);
            zos.putNextEntry(new ZipEntry(exportName + "-REMARKS.txt"));
            int nWarningCount = 0;
            BufferedReader in = new BufferedReader(new FileReader(warningFile));
            String sLine;
            while ((sLine = in.readLine()) != null) {
                zos.write((sLine + "\n").getBytes());
                progress.setCurrentStepProgress(nWarningCount++);
            }
            LOG.info("Number of Warnings for export (" + exportName + "): " + nWarningCount);
            in.close();
            zos.closeEntry();
        }
        warningFile.delete();

        zos.finish();
        zos.close();
        progress.setPercentageEnabled(true);
        progress.setCurrentStepProgress((short) 100);
    }

    public TreeMap<Integer, Comparable> writeGenotypeFile(OutputStream os, String sModule, Integer nAssemblyId, int nQueryChunkSize, MongoCollection<Document> varColl, Document varQuery, Map<String, String> markerSynonyms, File[] individualExportFiles, FileWriter warningFileWriter, ProgressIndicator progress) throws IOException, InterruptedException {
   		os.write(("# fjFile = GENOTYPE" + LINE_SEPARATOR).getBytes());
        
   		boolean fWorkingOnRuns = varColl.getNamespace().getCollectionName().equals(MongoTemplateManager.getMongoCollectionName(VariantRunData.class));
   		try (MongoCursor<Document> markerCursor = IExportHandler.getMarkerCursorWithCorrectCollation(varColl, varQuery, nAssemblyId, nQueryChunkSize)) {
	        while (markerCursor.hasNext()) {
	            Document exportVariant = markerCursor.next();
	            String markerId = (String) (fWorkingOnRuns ? Helper.readPossiblyNestedField(exportVariant, "_id." + VariantRunDataId.FIELDNAME_VARIANT_ID, ";", null) : exportVariant.get("_id"));
	            String exportedId = markerSynonyms == null ? markerId : markerSynonyms.get(markerId);
	            os.write(("\t" + exportedId).getBytes());
	        }
		}
		catch (IOException ioe) {
			LOG.error("Error occurred while writing Flapjack genotype file. Deleting " + individualExportFiles.length + " temporary individual-oriented files.", ioe);
			for (File f : individualExportFiles)
				f.delete();
		}
        os.write(LINE_SEPARATOR.getBytes());

        TreeMap<Integer, Comparable> problematicMarkerIndexToNameMap = new TreeMap<Integer, Comparable>();
        short nProgress = 0, nPreviousProgress = 0;

        int i = 0, nNConcurrentThreads = Math.max(1, Runtime.getRuntime().availableProcessors());	// use multiple threads so we can prepare several lines at once
        HashMap<Integer, StringBuffer> individualLines = new HashMap<>(nNConcurrentThreads);
        for (int j=0; j<nNConcurrentThreads; j++)
        	individualLines.put(j, new StringBuffer());
        final ArrayList<Thread> threadsToWaitFor = new ArrayList<>(nNConcurrentThreads);

        try
        {
	        for (final File f : individualExportFiles) {
	            if (progress.isAborted() || progress.getError() != null)
	            	return null;

	            final int nThreadIndex = i % nNConcurrentThreads;
	            Thread thread = new Thread() {
	                @Override
	                public void run() {
	    	        	StringBuffer indLine = individualLines.get(nThreadIndex);
	                	BufferedReader in = null;
	    	            try {
		    	            in = new BufferedReader(new FileReader(f));
	    	                String individualId, line = in.readLine();	// read sample id
	    	                if (line != null) {
	    	                    individualId = line;
	    	                    indLine.append(individualId);
	    	                } else {
	    	                    throw new Exception("Unable to read first line of temp export file " + f.getName());
	    	                }

	    	                int nMarkerIndex = 0;
	    	                while ((line = in.readLine()) != null) {
	    	                    List<String> genotypes = Helper.split(line, "|");
	    	                    HashMap<Object, Integer> genotypeCounts = new HashMap<Object, Integer>();	// will help us to keep track of missing genotypes
	    	                    int highestGenotypeCount = 0;
	    	                    String mostFrequentGenotype = null;
	    	                    for (String genotype : genotypes) {
	    	                        if (genotype == null) {
	    	                            continue;	/* skip missing genotypes */
	    	                        }
	    	
	    	                        int gtCount = 1 + Helper.getCountForKey(genotypeCounts, genotype);
	    	                        if (gtCount > highestGenotypeCount) {
	    	                            highestGenotypeCount = gtCount;
	    	                            mostFrequentGenotype = genotype;
	    	                        }
	    	                        genotypeCounts.put(genotype, gtCount);
	    	                    }
	    	
	    	                    if (genotypeCounts.size() > 1) {
	    	                    	if (warningFileWriter != null)
	    	                    		warningFileWriter.write("- Dissimilar genotypes found for variant n. " + nMarkerIndex + ", individual " + individualId + ". Exporting most frequent: " + mostFrequentGenotype + "\n");
	    	                        problematicMarkerIndexToNameMap.put(nMarkerIndex, "");
	    	                    }
	    	
	    	                    String[] alleles = mostFrequentGenotype == null ? new String[0] : mostFrequentGenotype.split(" ");
	    	                    if (alleles.length > 2) {
	    	                    	if (warningFileWriter != null)
	    	                    		warningFileWriter.write("- More than 2 alleles found for variant n. " + nMarkerIndex + ", individual " + individualId + ". Exporting only the first 2 alleles.\n");
	    	                        problematicMarkerIndexToNameMap.put(nMarkerIndex, "");
	    	                    }

	    	                    if (alleles.length == 0 || (alleles.length == 1 && alleles[0].length() == 0))
	    	                    	indLine.append("\t-");
	    	                    else
	    	                    {
	    		                    String all1 = alleles[0];
	    		                    String all2 = alleles[alleles.length == 1 ? 0 : 1];
	    		                    indLine.append("\t" + all1 + (!all2.equals(all1) ? "/" + all2 : ""));
	    	                    }

	    	                    nMarkerIndex++;
	    	                }
	    	            } catch (Exception e) {
	    	                LOG.error("Error exporting data", e);
	    	                progress.setError("Error exporting data: " + e.getClass().getSimpleName() + (e.getMessage() != null ? " - " + e.getMessage() : ""));
	    	            } finally {
	    	            	if (in != null)
								try {
									in.close();
								} catch (IOException e) {
									LOG.warn(e);
								}               
	    	            }

	    	            indLine.append(LINE_SEPARATOR);
	                }
	            };
	            
	            thread.start();
	            threadsToWaitFor.add(thread);

	            if (++i % nNConcurrentThreads == 0 || i == individualExportFiles.length) {
	                for (Thread t : threadsToWaitFor) // wait for all previously launched async threads
	               		t.join();
	                
	                for (int j=0; j<nNConcurrentThreads; j++) {
	                	os.write(individualLines.get(j).toString().getBytes());
	                	individualLines.put(j, new StringBuffer());
	                }

    	            nProgress = (short) (i * 100 / individualExportFiles.length);
    	            if (nProgress > nPreviousProgress) {
    	                progress.setCurrentStepProgress(nProgress);
    	                nPreviousProgress = nProgress;
    	            }
	                threadsToWaitFor.clear();
	            }
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
        return Arrays.asList(new String[]{"Exporting data to Flapjack format"});
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
    
	@Override
	public String[] getExportDataFileExtensions() {
		return new String[] {"genotype", "map"};
	}

	@Override
	public String getExportContentType() {
		return "application/x-fjzip";
	}

    @Override
    public int[] getSupportedPloidyLevels() {
        return null;
    }
}