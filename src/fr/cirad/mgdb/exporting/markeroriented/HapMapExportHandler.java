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

import fr.cirad.mgdb.exporting.tools.AsyncExportTool;
import fr.cirad.mgdb.exporting.tools.AsyncExportTool.AbstractDataOutputHandler;
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
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.springframework.data.mongodb.core.MongoTemplate;

import com.mongodb.DBCursor;
import com.mongodb.DBObject;

/**
 * The Class HapMapExportHandler.
 */
public class HapMapExportHandler extends AbstractMarkerOrientedExportHandler {

    /**
     * The Constant LOG.
     */
    private static final Logger LOG = Logger.getLogger(HapMapExportHandler.class);
    
    public static final String missingGenotype = "\tNN";

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
	
////    @Override
//    public void exportData_sync(OutputStream outputStream, String sModule, List<SampleId> sampleIDs1, List<SampleId> sampleIDs2, ProgressIndicator progress, DBCursor markerCursor, Map<Comparable, Comparable> markerSynonyms, HashMap<String, Integer> annotationFieldThresholds, HashMap<String, Integer> annotationFieldThresholds2, Map<SampleId, String> sampleIndexToIndividualMapToExport, Map<String, InputStream> readyToExportFiles) throws Exception {
//        byte[] missingGenotypeBytes = missingGenotype.getBytes();
//
//		List<String> individuals1 = getIndividualsFromSamples(sModule, sampleIDs1).stream().map(ind -> ind.getId()).collect(Collectors.toList());	
//		List<String> individuals2 = getIndividualsFromSamples(sModule, sampleIDs2).stream().map(ind -> ind.getId()).collect(Collectors.toList());
//
//		List<String> sortedIndividuals = sampleIndexToIndividualMapToExport.values().stream().distinct().sorted(new AlphaNumericComparator<String>()).collect(Collectors.toList());
//		
//		MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
//        File warningFile = File.createTempFile("export_warnings_", "");
//        FileWriter warningFileWriter = new FileWriter(warningFile);
//
//        int markerCount = markerCursor.count();
//
//        ZipOutputStream zos = new ZipOutputStream(outputStream);
//
//        if (readyToExportFiles != null) {
//            for (String readyToExportFile : readyToExportFiles.keySet()) {
//                zos.putNextEntry(new ZipEntry(readyToExportFile));
//                InputStream inputStream = readyToExportFiles.get(readyToExportFile);
//                byte[] dataBlock = new byte[1024];
//                int count = inputStream.read(dataBlock, 0, 1024);
//                while (count != -1) {
//                    zos.write(dataBlock, 0, count);
//                    count = inputStream.read(dataBlock, 0, 1024);
//                }
//                zos.closeEntry();
//            }
//        }
//
//        String exportName = sModule + "_" + markerCount + "variants_" + sortedIndividuals.size() + "individuals";
//        zos.putNextEntry(new ZipEntry(exportName + ".hapmap"));
//        String header = "rs#" + "\t" + "alleles" + "\t" + "chrom" + "\t" + "pos" + "\t" + "strand" + "\t" + "assembly#" + "\t" + "center" + "\t" + "protLSID" + "\t" + "assayLSID" + "\t" + "panelLSID" + "\t" + "QCcode";
//        zos.write(header.getBytes());
//        for (String individual : sortedIndividuals)
//            zos.write(("\t" + individual).getBytes());
//        zos.write((LINE_SEPARATOR).getBytes());
//
//        int avgObjSize = (Integer) mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantRunData.class)).getStats().get("avgObjSize");
//        int nChunkSize = nMaxChunkSizeInMb * 1024 * 1024 / avgObjSize;
//        short nProgress = 0, nPreviousProgress = 0;
//        long nLoadedMarkerCount = 0;
//
//        while (markerCursor == null || markerCursor.hasNext()) {
//            int nLoadedMarkerCountInLoop = 0;
//            Map<Comparable, String> markerChromosomalPositions = new LinkedHashMap<Comparable, String>();
//            boolean fStartingNewChunk = true;
//            markerCursor.batchSize(nChunkSize);
//            while (markerCursor.hasNext() && (fStartingNewChunk || nLoadedMarkerCountInLoop % nChunkSize != 0)) {
//                DBObject exportVariant = markerCursor.next();
//                DBObject refPos = (DBObject) exportVariant.get(VariantData.FIELDNAME_REFERENCE_POSITION);
//                markerChromosomalPositions.put((Comparable) exportVariant.get("_id"), refPos == null ? null : (refPos.get(ReferencePosition.FIELDNAME_SEQUENCE) + ":" + refPos.get(ReferencePosition.FIELDNAME_START_SITE)));
//                nLoadedMarkerCountInLoop++;
//                fStartingNewChunk = false;
//            }
//
//            List<Comparable> currentMarkers = new ArrayList<Comparable>(markerChromosomalPositions.keySet());
//            LinkedHashMap<VariantData, Collection<VariantRunData>> variantsAndRuns = MgdbDao.getSampleGenotypes(mongoTemplate, sampleIndexToIndividualMapToExport.keySet(), currentMarkers, true, null /*new Sort(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ChromosomalPosition.FIELDNAME_SEQUENCE).and(new Sort(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ChromosomalPosition.FIELDNAME_START_SITE))*/);	// query mongo db for matching genotypes			
//            for (VariantData variant : variantsAndRuns.keySet()) // read data and write results into temporary files (one per sample)
//            {
//                Comparable variantId = variant.getId();
//                if (markerSynonyms != null) {
//                    Comparable syn = markerSynonyms.get(variantId);
//                    if (syn != null) {
//                        variantId = syn;
//                    }
//                }
//
//                boolean fIsSNP = variant.getType().equals(Type.SNP.toString());
//
//                String refPos = markerChromosomalPositions.get(variant.getId());
//                String[] chromAndPos = refPos == null ? null : refPos.split(":");
//                zos.write(((variantId == null ? variant.getId() : variantId) + "\t" + StringUtils.join(variant.getKnownAlleleList(), "/") + "\t" + (chromAndPos == null ? 0 : chromAndPos[0]) + "\t" + (chromAndPos == null ? 0 : Long.parseLong(chromAndPos[1])) + "\t" + "+").getBytes());
//                for (int j = 0; j < 6; j++) {
//                    zos.write(("\t" + "NA").getBytes());
//                }
//
//                Map<String, List<String>> individualGenotypes = new TreeMap<String, List<String>>(new AlphaNumericComparator<String>());
//                Collection<VariantRunData> runs = variantsAndRuns.get(variant);
//                if (runs != null) {
//                    for (VariantRunData run : runs) {
//                        for (Integer sampleIndex : run.getSampleGenotypes().keySet()) {
//                            SampleGenotype sampleGenotype = run.getSampleGenotypes().get(sampleIndex);
//                            String individualId = sampleIndexToIndividualMapToExport.get(new SampleId(run.getId().getProjectId(), sampleIndex));
//                            
//							if (!VariantData.gtPassesVcfAnnotationFilters(individualId, sampleGenotype, individuals1, annotationFieldThresholds, individuals2, annotationFieldThresholds2))
//								continue;	// skip genotype
//							
//                            String gtCode = run.getSampleGenotypes().get(sampleIndex).getCode();
//                            List<String> storedIndividualGenotypes = individualGenotypes.get(individualId);
//                            if (storedIndividualGenotypes == null) {
//                                storedIndividualGenotypes = new ArrayList<String>();
//                                individualGenotypes.put(individualId, storedIndividualGenotypes);
//                            }
//                            storedIndividualGenotypes.add(gtCode);
//                        }
//                    }
//                }
//
//                int writtenGenotypeCount = 0;
//                for (String individual : individualGenotypes.keySet() /* we use this list because it has the proper ordering */) {
//                    int individualIndex = sortedIndividuals.indexOf(individual);
//                    while (writtenGenotypeCount < individualIndex) {
//                        zos.write(missingGenotypeBytes);
//                        writtenGenotypeCount++;
//                    }
//
//                    List<String> genotypes = individualGenotypes.get(individual);
//                    HashMap<Object, Integer> genotypeCounts = new HashMap<Object, Integer>();	// will help us to keep track of missing genotypes
//                    int highestGenotypeCount = 0;
//                    String mostFrequentGenotype = null;
//                    if (genotypes != null) {
//                        for (String genotype : genotypes) {
//                            if (genotype.length() == 0) {
//                                continue;	/* skip missing genotypes */
//                            }
//
//                            int gtCount = 1 + Helper.getCountForKey(genotypeCounts, genotype);
//                            if (gtCount > highestGenotypeCount) {
//                                highestGenotypeCount = gtCount;
//                                mostFrequentGenotype = genotype;
//                            }
//                            genotypeCounts.put(genotype, gtCount);
//                        }
//                    }
//
//                    byte[] exportedGT = mostFrequentGenotype == null ? missingGenotypeBytes : ("\t" + StringUtils.join(variant.getAllelesFromGenotypeCode(mostFrequentGenotype), fIsSNP ? "" : "/")).getBytes();
//                    zos.write(exportedGT);
//                    writtenGenotypeCount++;
//
//                    if (genotypeCounts.size() > 1) {
//                        warningFileWriter.write("- Dissimilar genotypes found for variant " + (variantId == null ? variant.getId() : variantId) + ", individual " + individual + ". Exporting most frequent: " + new String(exportedGT) + "\n");
//                    }
//                }
//
//                while (writtenGenotypeCount < sortedIndividuals.size()) {
//                    zos.write(missingGenotypeBytes);
//                    writtenGenotypeCount++;
//                }
//                zos.write((LINE_SEPARATOR).getBytes());
//            }
//
//            if (progress.hasAborted()) {
//                warningFileWriter.close();
//                return;
//            }
//
//            nLoadedMarkerCount += nLoadedMarkerCountInLoop;
//            nProgress = (short) (nLoadedMarkerCount * 100 / markerCount);
//            if (nProgress > nPreviousProgress) {
//                //				if (nProgress%5 == 0)
//                //					LOG.info("========================= exportData: " + nProgress + "% =========================" + (System.currentTimeMillis() - before)/1000 + "s");
//                progress.setCurrentStepProgress(nProgress);
//                nPreviousProgress = nProgress;
//            }
//        }
//        zos.closeEntry();
//
//		warningFileWriter.close();
//		if (warningFile.length() > 0) {
//			zos.putNextEntry(new ZipEntry(exportName + "-REMARKS.txt"));
//			int nWarningCount = 0;
//			BufferedReader in = new BufferedReader(new FileReader(warningFile));
//			String sLine;
//			while ((sLine = in.readLine()) != null) {
//				zos.write((sLine + "\n").getBytes());
//				nWarningCount++;
//			}
//			LOG.info("Number of Warnings for export (" + exportName + "): " + nWarningCount);
//			in.close();
//			zos.closeEntry();
//		}
//        warningFile.delete();
//
//        zos.finish();
//        zos.close();
//        progress.setCurrentStepProgress((short) 100);
//    }

    /* (non-Javadoc)
 * @see fr.cirad.mgdb.exporting.markeroriented.AbstractMarkerOrientedExportHandler#exportData(java.io.OutputStream, java.lang.String, java.util.List, fr.cirad.tools.ProgressIndicator, com.mongodb.DBCursor, java.util.Map, int, int, java.util.Map)
     */
    @Override
    public void exportData(OutputStream outputStream, String sModule, List<SampleId> sampleIDs1, List<SampleId> sampleIDs2, ProgressIndicator progress, DBCursor markerCursor, Map<Comparable, Comparable> markerSynonyms, HashMap<String, Integer> annotationFieldThresholds, HashMap<String, Integer> annotationFieldThresholds2, Map<SampleId, String> sampleIndexToIndividualMapToExport, Map<String, InputStream> readyToExportFiles) throws Exception {
        
		List<String> individuals1 = getIndividualsFromSamples(sModule, sampleIDs1).stream().map(ind -> ind.getId()).collect(Collectors.toList());	
		List<String> individuals2 = getIndividualsFromSamples(sModule, sampleIDs2).stream().map(ind -> ind.getId()).collect(Collectors.toList());

		ArrayList<SampleId> sampleIDs = (ArrayList<SampleId>) CollectionUtils.union(sampleIDs1, sampleIDs2);
		Collection<Individual> individuals = getIndividualsFromSamples(sModule, sampleIDs);
		LinkedHashMap<SampleId, String> sampleIdToIndividualMap = new LinkedHashMap<SampleId, String>();
		for (int i=0; i<sampleIDs.size(); i++)
			sampleIdToIndividualMap.put(sampleIDs.get(i), ((List<Individual>)individuals).get(i).getId());
		
		// from now on it will have the proper ordering (until then we needed the ordering to remain consistent wtih that of SampleIDs)
		individuals = new TreeSet<Individual>(individuals);
		List<String> sortedIndividuals = individuals.stream().map(ind -> ind.getId()).distinct().collect(Collectors.toList());
		
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

        StringBuffer indSB = new StringBuffer();
        String exportName = sModule + "__" + markerCount + "variants__" + sortedIndividuals.size() + "individuals";
        zos.putNextEntry(new ZipEntry(exportName + ".hapmap"));
        String header = "rs#" + "\t" + "alleles" + "\t" + "chrom" + "\t" + "pos" + "\t" + "strand" + "\t" + "assembly#" + "\t" + "center" + "\t" + "protLSID" + "\t" + "assayLSID" + "\t" + "panelLSID" + "\t" + "QCcode";
        indSB.append(header);
        for (String individual : sortedIndividuals)
        	indSB.append("\t" + individual);
        indSB.append(LINE_SEPARATOR);

		AbstractDataOutputHandler<Integer, LinkedHashMap<VariantData, Collection<VariantRunData>>> dataOutputHandler = new AbstractDataOutputHandler<Integer, LinkedHashMap<VariantData, Collection<VariantRunData>>>() {				
			@Override
			public Void call() {
				StringBuffer sb = new StringBuffer();
				for (VariantData variant : variantDataChunkMap.keySet())
					try
					{
						Comparable variantId = variant.getId();
		                if (markerSynonyms != null) {
		                	Comparable syn = markerSynonyms.get(variantId);
		                    if (syn != null) {
		                        variantId = syn;
		                    }
		                }
	
		                boolean fIsSNP = variant.getType().equals(Type.SNP.toString());
	
						sb.append((variantId == null ? variant.getId() : variantId) + "\t" + StringUtils.join(variant.getKnownAlleleList(), "/") + "\t" + (variant.getReferencePosition() == null ? 0 : variant.getReferencePosition().getSequence()) + "\t" + (variant.getReferencePosition() == null ? 0 : variant.getReferencePosition().getStartSite()) + "\t" + "+\tNA\tNA\tNA\tNA\tNA\tNA");
	
		                Map<String, List<String>> individualGenotypes = new TreeMap<String, List<String>>(new AlphaNumericComparator<String>());
		                Collection<VariantRunData> runs = variantDataChunkMap.get(variant);
		                if (runs != null) {
		                    for (VariantRunData run : runs) {
								for (Integer sampleIndex : run.getSampleGenotypes().keySet()) {
									SampleGenotype sampleGenotype = run.getSampleGenotypes().get(sampleIndex);
									String individualId = sampleIdToIndividualMap.get(new SampleId(run.getId().getProjectId(), sampleIndex));
		                            
									if (!VariantData.gtPassesVcfAnnotationFilters(individualId, sampleGenotype, individuals1, annotationFieldThresholds, individuals2, annotationFieldThresholds2))
										continue;	// skip genotype
									
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
	
		                int writtenGenotypeCount = 0;
		                for (String individual : individualGenotypes.keySet() /* we use this list because it has the proper ordering */) {
		                    int individualIndex = sortedIndividuals.indexOf(individual);
		                    while (writtenGenotypeCount < individualIndex) {
		                        sb.append(missingGenotype);
		                        writtenGenotypeCount++;
		                    }
	
		                    List<String> genotypes = individualGenotypes.get(individual);
		                    HashMap<Object, Integer> genotypeCounts = new HashMap<Object, Integer>();	// will help us to keep track of missing genotypes
		                    int highestGenotypeCount = 0;
		                    String mostFrequentGenotype = null;
		                    if (genotypes != null) {
		                        for (String genotype : genotypes) {
		                            if (genotype == null)
		                                continue;	/* skip missing genotypes */
	
		                            int gtCount = 1 + Helper.getCountForKey(genotypeCounts, genotype);
		                            if (gtCount > highestGenotypeCount) {
		                                highestGenotypeCount = gtCount;
		                                mostFrequentGenotype = genotype;
		                            }
		                            genotypeCounts.put(genotype, gtCount);
		                        }
		                    }
	
		                    String exportedGT = mostFrequentGenotype == null || mostFrequentGenotype.isEmpty() ? missingGenotype : ("\t" + StringUtils.join(variant.getAllelesFromGenotypeCode(mostFrequentGenotype), fIsSNP ? "" : "/"));
		                    sb.append(exportedGT);
		                    writtenGenotypeCount++;
	
		                    if (genotypeCounts.size() > 1) {
		                        warningFileWriter.write("- Dissimilar genotypes found for variant " + (variantId == null ? variant.getId() : variantId) + ", individual " + individual + ". Exporting most frequent: " + new String(exportedGT) + "\n");
		                    }
		                }
	
		                while (writtenGenotypeCount < sortedIndividuals.size()) {
		                    sb.append(missingGenotype);
		                    writtenGenotypeCount++;
		                }
		                sb.append(LINE_SEPARATOR);
	
			            if (progress.hasAborted()) {
			                warningFileWriter.close();
			                return null;
			            }
	                }
					catch (Exception e)
					{
						if (progress.getError() == null)	// only log this once
							LOG.debug("Unable to export " + variant.getId(), e);
						progress.setError("Unable to export " + variant.getId() + ": " + e.getMessage());
					}
				
                try
                {
    				zos.write(sb.toString().getBytes());
				}
                catch (IOException ioe)
                {
                	progress.setError("Unable to export data for " + variantDataChunkMap.keySet() + ": " + ioe.getMessage());
                }
				return null;
			}
		};

    	Number avgObjSize = (Number) mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantRunData.class)).getStats().get("avgObjSize");
		int nQueryChunkSize = (int) Math.max(1, (nMaxChunkSizeInMb*1024*1024 / avgObjSize.doubleValue()) / AsyncExportTool.WRITING_QUEUE_CAPACITY);
		AsyncExportTool syncExportTool = new AsyncExportTool(markerCursor, markerCount, nQueryChunkSize, mongoTemplate, sampleIDs, dataOutputHandler, progress);
		syncExportTool.launch();

		while (progress.getCurrentStepProgress() < 100 && !progress.hasAborted())
			Thread.sleep(500);
		
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
