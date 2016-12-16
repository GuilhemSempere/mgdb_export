/*******************************************************************************
 * MGDB Export - Mongo Genotype DataBase, export handlers
 * Copyright (C) 2016 <CIRAD>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License, version 3 as
 * published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * See <http://www.gnu.org/licenses/agpl.html> for details about
 * GNU Affero General Public License V3.
 *******************************************************************************/
package fr.cirad.mgdb.exporting.markeroriented;

import fr.cirad.mgdb.model.mongo.maintypes.Sequence;
import fr.cirad.mgdb.model.mongo.maintypes.SequenceStats;
import fr.cirad.mgdb.exporting.tools.AsyncExportTool;
import fr.cirad.mgdb.exporting.tools.AsyncExportTool.AbstractDataOutputHandler;
import fr.cirad.mgdb.model.mongo.maintypes.DBVCFHeader;
import fr.cirad.mgdb.model.mongo.maintypes.DBVCFHeader.VcfHeaderId;
import fr.cirad.mgdb.model.mongo.maintypes.Individual;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongo.subtypes.SampleId;
import fr.cirad.mgdb.model.mongodao.MgdbDao;
import fr.cirad.tools.AlphaNumericStringComparator;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mongo.MongoTemplateManager;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.VariantContext;
//import htsjdk.variant.variantcontext.writer.AsyncVariantContextWriter;
import htsjdk.variant.variantcontext.writer.CustomVCFWriter;
//import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
//import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

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
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.log4j.Logger;
import org.bson.types.ObjectId;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.mapreduce.MapReduceResults;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;

import com.mongodb.BasicDBObject;
import com.mongodb.DBCursor;
import com.mongodb.DBObject;

/**
 * The Class VcfExportHandler.
 */
public class VcfExportHandler extends AbstractMarkerOrientedExportHandler {

	/** The Constant LOG. */
	private static final Logger LOG = Logger.getLogger(VcfExportHandler.class);

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatName()
	 */
	@Override
	public String getExportFormatName()
	{
		return "VCF";
	}

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
	 */
	@Override
	public String getExportFormatDescription()
	{
		return "Exports data in Variant Call Format. See <a target='_blank' href='http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41'>http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41</a> for more details.";
	}

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
	 */
	@Override
	public List<String> getStepList()
	{
		return Arrays.asList(new String[] {"Creating sequence list", "Exporting data to VCF format"});
	}

	/**
	 * Creates the sam sequence dictionary.
	 *
	 * @param sModule the module
	 * @param sequenceIDs the sequence i ds
	 * @return the SAM sequence dictionary
	 * @throws Exception the exception
	 */
	public SAMSequenceDictionary createSAMSequenceDictionary(String sModule, Collection<String> sequenceIDs) throws Exception
	{
		SAMSequenceDictionary dict1 = new SAMSequenceDictionary();
		MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
		String sequenceSeqCollName = MongoTemplateManager.getMongoCollectionName(Sequence.class);
		if (mongoTemplate.collectionExists(sequenceSeqCollName) && sequenceIDs.size() > 1)
		{
			long before = System.currentTimeMillis();
			Query query = new Query(Criteria.where("_id").in(sequenceIDs));
			String mapFunction = "function() { emit(this._id, this." + Sequence.FIELDNAME_SEQUENCE + ".length); }";
			String reduceFunction = "function() { }";
			MapReduceResults<Map> rs = MongoTemplateManager.get(sModule).mapReduce(query, sequenceSeqCollName, mapFunction, reduceFunction, Map.class);
			Iterator<Map> it = rs.iterator();
			while (it.hasNext())
			{
				Map map = it.next();
				dict1.addSequence(new SAMSequenceRecord((String) map.get("_id"), ((Double) map.get("value")).intValue()));
			}
	    	LOG.info("createSAMSequenceDictionary took " + (System.currentTimeMillis() - before)/1000d + "s to write " + sequenceIDs.size() + " sequences");
		}
		else
			LOG.info("No sequence data was found. No SAMSequenceDictionary could be created.");
	    return dict1;
	}

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.markeroriented.AbstractMarkerOrientedExportHandler#exportData(java.io.OutputStream, java.lang.String, java.util.List, fr.cirad.tools.ProgressIndicator, com.mongodb.DBCursor, java.util.Map, int, int, java.util.Map)
	 */
	@Override
	public void exportData(OutputStream outputStream, String sModule, List<SampleId> sampleIDs, ProgressIndicator progress, DBCursor markerCursor, Map<Comparable, Comparable> markerSynonyms, int nMinimumGenotypeQuality, int nMinimumReadDepth, Map<String, InputStream> readyToExportFiles) throws Exception
	{
		Integer projectId = null;
		for (SampleId spId : sampleIDs)
		{
			if (projectId == null)
				projectId = spId.getProject();
			else if (projectId != spId.getProject())
			{
				projectId = 0;
				break;	// more than one project are involved: no header will be written
			}
		}

		File warningFile = File.createTempFile("export_warnings_", "");
		FileWriter warningFileWriter = new FileWriter(warningFile);

		MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
		int markerCount = markerCursor.count();
		ZipOutputStream zos = new ZipOutputStream(outputStream);

		if (readyToExportFiles != null)
			for (String readyToExportFile : readyToExportFiles.keySet())
			{
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

		LinkedHashMap<SampleId, String> sampleIDToIndividualIdMap = new LinkedHashMap<SampleId, String>();
		ArrayList<String> individualList = new ArrayList<String>();
		List<Individual> individuals = getIndividualsFromSamples(sModule, sampleIDs);
		for (int i=0; i<sampleIDs.size(); i++) {
			String individualId = individuals.get(i).getId();
			sampleIDToIndividualIdMap.put(sampleIDs.get(i), individualId);
			if (!individualList.contains(individualId)) {
				individualList.add(individualId);
			}
		}

		String exportName = sModule + "_" + markerCount + "variants_" + individualList.size() + "individuals";

		int avgObjSize = (Integer) mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantRunData.class)).getStats().get("avgObjSize");
		int nQueryChunkSize = nMaxChunkSizeInMb*1024*1024 / avgObjSize;

		VariantContextWriter writer = null;
		try
		{
			zos.putNextEntry(new ZipEntry(exportName + ".vcf"));
			List<String> distinctSequenceNames = new ArrayList<String>();

			String sequenceSeqCollName = MongoTemplateManager.getMongoCollectionName(Sequence.class);
			if (mongoTemplate.collectionExists(sequenceSeqCollName))
			{
				DBCursor markerCursorCopy = markerCursor.copy();
				markerCursorCopy.batchSize(nQueryChunkSize);
				while (markerCursorCopy.hasNext())
				{
					int nLoadedMarkerCountInLoop = 0;
					boolean fStartingNewChunk = true;
					while (markerCursorCopy.hasNext() && (fStartingNewChunk || nLoadedMarkerCountInLoop%nQueryChunkSize != 0)) {
						DBObject exportVariant = markerCursorCopy.next();
						String chr = (String) ((DBObject) exportVariant.get(VariantData.FIELDNAME_REFERENCE_POSITION)).get(ReferencePosition.FIELDNAME_SEQUENCE);
						if (!distinctSequenceNames.contains(chr))
							distinctSequenceNames.add(chr);
					}
				}
				markerCursorCopy.close();
			}

			Collections.sort(distinctSequenceNames, new AlphaNumericStringComparator());
			SAMSequenceDictionary dict = createSAMSequenceDictionary(sModule, distinctSequenceNames);
			writer = new CustomVCFWriter(null, zos, dict, false, false, true);
//			VariantContextWriterBuilder vcwb = new VariantContextWriterBuilder();
//			vcwb.unsetOption(Options.INDEX_ON_THE_FLY);
//			vcwb.unsetOption(Options.DO_NOT_WRITE_GENOTYPES);
//			vcwb.setOption(Options.USE_ASYNC_IOINDEX_ON_THE_FLY);
//			vcwb.setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
//			vcwb.setReferenceDictionary(dict);
//			writer = vcwb.build();
//			writer = new AsyncVariantContextWriter(writer, 3000);

			progress.moveToNextStep();	// done with dictionary
			DBCursor headerCursor = mongoTemplate.getCollection(MongoTemplateManager.getMongoCollectionName(DBVCFHeader.class)).find(new BasicDBObject("_id." + VcfHeaderId.FIELDNAME_PROJECT, projectId));
			Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
			boolean fWriteCommandLine = true, fWriteEngineHeaders = true;	// default values

			while (headerCursor.hasNext())
			{
				DBVCFHeader dbVcfHeader = DBVCFHeader.fromDBObject(headerCursor.next());
				headerLines.addAll(dbVcfHeader.getHeaderLines());

				// Add sequence header lines (not stored in our vcf header collection)
				BasicDBObject projection = new BasicDBObject(SequenceStats.FIELDNAME_SEQUENCE_LENGTH, true);
				int nSequenceIndex = 0;
				for (String sequenceName : distinctSequenceNames)
				{
					String sequenceInfoCollName = MongoTemplateManager.getMongoCollectionName(SequenceStats.class);
					boolean fCollectionExists = mongoTemplate.collectionExists(sequenceInfoCollName);
					if (fCollectionExists) {
						DBObject record = mongoTemplate.getCollection(sequenceInfoCollName).findOne(new Query(Criteria.where("_id").is(sequenceName)).getQueryObject(), projection);
						if (record == null)
						{
							LOG.warn("Sequence '" + sequenceName + "' not found in collection " + sequenceInfoCollName);
							continue;
						}

						Map<String, String> sequenceLineData = new LinkedHashMap<String, String>();
						sequenceLineData.put("ID", (String) record.get("_id"));
						sequenceLineData.put("length", ((Number) record.get(SequenceStats.FIELDNAME_SEQUENCE_LENGTH)).toString());
						headerLines.add(new VCFContigHeaderLine(sequenceLineData, nSequenceIndex++));
					}
				}
				fWriteCommandLine = headerCursor.size() == 1 && dbVcfHeader.getWriteCommandLine();	// wouldn't make sense to include command lines for several runs
				if (!dbVcfHeader.getWriteEngineHeaders())
					fWriteEngineHeaders = false;
			}
			headerCursor.close();

			VCFHeader header = new VCFHeader(headerLines, individualList);
			header.setWriteCommandLine(fWriteCommandLine);
			header.setWriteEngineHeaders(fWriteEngineHeaders);
			writer.writeHeader(header);

			short nProgress = 0, nPreviousProgress = 0;
			long nLoadedMarkerCount = 0;
			HashMap<SampleId, Comparable /*phID*/> phasingIDsBySample = new HashMap<SampleId, Comparable>();

			while (markerCursor.hasNext())
			{
				if (progress.hasAborted())
					return;

				int nLoadedMarkerCountInLoop = 0;
				boolean fStartingNewChunk = true;
				markerCursor.batchSize(nQueryChunkSize);
				List<Comparable> currentMarkers = new ArrayList<Comparable>();
				while (markerCursor.hasNext() && (fStartingNewChunk || nLoadedMarkerCountInLoop%nQueryChunkSize != 0))
				{
					DBObject exportVariant = markerCursor.next();
					currentMarkers.add((Comparable) exportVariant.get("_id"));
					nLoadedMarkerCountInLoop++;
					fStartingNewChunk = false;
				}

				LinkedHashMap<VariantData, Collection<VariantRunData>> variantsAndRuns = MgdbDao.getSampleGenotypes(mongoTemplate, sampleIDs, currentMarkers, true, null /*new Sort(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ChromosomalPosition.FIELDNAME_SEQUENCE).and(new Sort(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ChromosomalPosition.FIELDNAME_START_SITE))*/);	// query mongo db for matching genotypes
				for (VariantData variant : variantsAndRuns.keySet())
				{
					VariantContext vc = variant.toVariantContext(variantsAndRuns.get(variant), !ObjectId.isValid(variant.getId().toString()), sampleIDToIndividualIdMap, phasingIDsBySample, nMinimumGenotypeQuality, nMinimumReadDepth, warningFileWriter, markerSynonyms == null ? variant.getId() : markerSynonyms.get(variant.getId()));
					try
					{
						writer.add(vc);
					}
					catch (Throwable t)
					{
						Exception e = new Exception("Unable to convert to VariantContext: " + variant.getId(), t);
						LOG.debug("error", e);
						throw e;
					}

					if (nLoadedMarkerCountInLoop > currentMarkers.size())
						LOG.error("Bug: writing variant number " + nLoadedMarkerCountInLoop + " (only " + currentMarkers.size() + " variants expected)");
				}

	            nLoadedMarkerCount += nLoadedMarkerCountInLoop;
				nProgress = (short) (nLoadedMarkerCount * 100 / markerCount);
				if (nProgress > nPreviousProgress) {
					progress.setCurrentStepProgress(nProgress);
					nPreviousProgress = nProgress;
				}
			}
			progress.setCurrentStepProgress((short) 100);
			zos.closeEntry();
		}
		catch (Exception e)
		{
			LOG.error("Error exporting", e);
			progress.setError(e.getMessage());
			return;
		}
		finally
		{
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
			if (writer != null)
				try
				{
					zos.finish();
					writer.close();
				}
				catch (Throwable ignored)
				{}
		}
	}
	
	public void exportData_new(OutputStream outputStream, String sModule, List<SampleId> sampleIDs, ProgressIndicator progress, DBCursor markerCursor, Map<Comparable, Comparable> markerSynonyms, int nMinimumGenotypeQuality, int nMinimumReadDepth, Map<String, InputStream> readyToExportFiles) throws Exception
	{
		Integer projectId = null;
		for (SampleId spId : sampleIDs)
		{
			if (projectId == null)
				projectId = spId.getProject();
			else if (projectId != spId.getProject())
			{
				projectId = 0;
				break;	// more than one project are involved: no header will be written
			}
		}

		File warningFile = File.createTempFile("export_warnings_", "");
		FileWriter warningFileWriter = new FileWriter(warningFile);

		MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
		int markerCount = markerCursor.count();
		ZipOutputStream zos = new ZipOutputStream(outputStream);

		if (readyToExportFiles != null)
			for (String readyToExportFile : readyToExportFiles.keySet())
			{
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

		LinkedHashMap<SampleId, String> sampleIDToIndividualIdMap = new LinkedHashMap<SampleId, String>();
		ArrayList<String> individualList = new ArrayList<String>();
		List<Individual> individuals = getIndividualsFromSamples(sModule, sampleIDs);
		for (int i=0; i<sampleIDs.size(); i++) {
			String individualId = individuals.get(i).getId();
			sampleIDToIndividualIdMap.put(sampleIDs.get(i), individualId);
			if (!individualList.contains(individualId)) {
				individualList.add(individualId);
			}
		}

		String exportName = sModule + "_" + markerCount + "variants_" + individualList.size() + "individuals";

		int avgObjSize = (Integer) mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantRunData.class)).getStats().get("avgObjSize");
		int nQueryChunkSize = nMaxChunkSizeInMb*1024*1024 / avgObjSize;

		VariantContextWriter writer = null;
		try
		{
			zos.putNextEntry(new ZipEntry(exportName + ".vcf"));
			List<String> distinctSequenceNames = new ArrayList<String>();

			String sequenceSeqCollName = MongoTemplateManager.getMongoCollectionName(Sequence.class);
			if (mongoTemplate.collectionExists(sequenceSeqCollName))
			{
				DBCursor markerCursorCopy = markerCursor.copy();
				markerCursorCopy.batchSize(nQueryChunkSize);
				while (markerCursorCopy.hasNext())
				{
					int nLoadedMarkerCountInLoop = 0;
					boolean fStartingNewChunk = true;
					while (markerCursorCopy.hasNext() && (fStartingNewChunk || nLoadedMarkerCountInLoop%nQueryChunkSize != 0)) {
						DBObject exportVariant = markerCursorCopy.next();
						String chr = (String) ((DBObject) exportVariant.get(VariantData.FIELDNAME_REFERENCE_POSITION)).get(ReferencePosition.FIELDNAME_SEQUENCE);
						if (!distinctSequenceNames.contains(chr))
							distinctSequenceNames.add(chr);
					}
				}
				markerCursorCopy.close();
			}

			Collections.sort(distinctSequenceNames, new AlphaNumericStringComparator());
			SAMSequenceDictionary dict = createSAMSequenceDictionary(sModule, distinctSequenceNames);
			writer = new CustomVCFWriter(null, zos, dict, false, false, true);
//			VariantContextWriterBuilder vcwb = new VariantContextWriterBuilder();
//			vcwb.unsetOption(Options.INDEX_ON_THE_FLY);
//			vcwb.unsetOption(Options.DO_NOT_WRITE_GENOTYPES);
//			vcwb.setOption(Options.USE_ASYNC_IOINDEX_ON_THE_FLY);
//			vcwb.setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
//			vcwb.setReferenceDictionary(dict);
//			writer = vcwb.build();
//			writer = new AsyncVariantContextWriter(writer, 3000);

			progress.moveToNextStep();	// done with dictionary
			DBCursor headerCursor = mongoTemplate.getCollection(MongoTemplateManager.getMongoCollectionName(DBVCFHeader.class)).find(new BasicDBObject("_id." + VcfHeaderId.FIELDNAME_PROJECT, projectId));
			Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
			boolean fWriteCommandLine = true, fWriteEngineHeaders = true;	// default values

			while (headerCursor.hasNext())
			{
				DBVCFHeader dbVcfHeader = DBVCFHeader.fromDBObject(headerCursor.next());
				headerLines.addAll(dbVcfHeader.getHeaderLines());

				// Add sequence header lines (not stored in our vcf header collection)
				BasicDBObject projection = new BasicDBObject(SequenceStats.FIELDNAME_SEQUENCE_LENGTH, true);
				int nSequenceIndex = 0;
				for (String sequenceName : distinctSequenceNames)
				{
					String sequenceInfoCollName = MongoTemplateManager.getMongoCollectionName(SequenceStats.class);
					boolean fCollectionExists = mongoTemplate.collectionExists(sequenceInfoCollName);
					if (fCollectionExists) {
						DBObject record = mongoTemplate.getCollection(sequenceInfoCollName).findOne(new Query(Criteria.where("_id").is(sequenceName)).getQueryObject(), projection);
						if (record == null)
						{
							LOG.warn("Sequence '" + sequenceName + "' not found in collection " + sequenceInfoCollName);
							continue;
						}

						Map<String, String> sequenceLineData = new LinkedHashMap<String, String>();
						sequenceLineData.put("ID", (String) record.get("_id"));
						sequenceLineData.put("length", ((Number) record.get(SequenceStats.FIELDNAME_SEQUENCE_LENGTH)).toString());
						headerLines.add(new VCFContigHeaderLine(sequenceLineData, nSequenceIndex++));
					}
				}
				fWriteCommandLine = headerCursor.size() == 1 && dbVcfHeader.getWriteCommandLine();	// wouldn't make sense to include command lines for several runs
				if (!dbVcfHeader.getWriteEngineHeaders())
					fWriteEngineHeaders = false;
			}
			headerCursor.close();

			VCFHeader header = new VCFHeader(headerLines, individualList);
			header.setWriteCommandLine(fWriteCommandLine);
			header.setWriteEngineHeaders(fWriteEngineHeaders);
			writer.writeHeader(header);

			HashMap<SampleId, Comparable /*phID*/> phasingIDsBySample = new HashMap<SampleId, Comparable>();
			
			final VariantContextWriter finalVariantContextWriter = writer;
			AbstractDataOutputHandler<Integer, LinkedHashMap<VariantData, Collection<VariantRunData>>[]> dataOutputHandler = new AbstractDataOutputHandler<Integer, LinkedHashMap<VariantData, Collection<VariantRunData>>[]>() {				
				@Override
				public Void call() {
					LinkedHashMap<VariantData, Collection<VariantRunData>> variantsAndRuns = variantDataArray[chunkIndex];
					for (VariantData variant : variantsAndRuns.keySet())
					{
						if (!progress.hasAborted())
							try
							{
								VariantContext vc = variant.toVariantContext(variantsAndRuns.get(variant), !ObjectId.isValid(variant.getId().toString()), sampleIDToIndividualIdMap, phasingIDsBySample, nMinimumGenotypeQuality, nMinimumReadDepth, warningFileWriter, markerSynonyms == null ? variant.getId() : markerSynonyms.get(variant.getId()));
								finalVariantContextWriter.add(vc);
							}
							catch (Exception e)
							{
								LOG.debug("Unable to convert to VariantContext", e);
								progress.setError("Unable to convert to VariantContext: " + e);
							}
					}
					return null;
				}
			};			
			
			AsyncExportTool syncExportTool = new AsyncExportTool(markerCursor, markerCount, nQueryChunkSize, mongoTemplate, sampleIDs, dataOutputHandler, progress);
			syncExportTool.launch();

			while (progress.getCurrentStepProgress() < 100 && !progress.hasAborted())
				Thread.sleep(500);
			zos.closeEntry();
		}
		catch (Exception e)
		{
			LOG.error("Error exporting", e);
			progress.setError(e.getMessage());
			return;
		}
		finally
		{
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
			if (writer != null)
				try
				{
					zos.finish();
					writer.close();
				}
				catch (Throwable ignored)
				{}
		}
	}
}