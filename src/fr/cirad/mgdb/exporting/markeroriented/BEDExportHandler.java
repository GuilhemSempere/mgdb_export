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
 * See <http://www.gnu.org/licenses/gpl-3.0.html> for details about
 * GNU Affero General Public License V3.
 *******************************************************************************/
package fr.cirad.mgdb.exporting.markeroriented;

import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.log4j.Logger;
import org.springframework.data.mongodb.core.MongoTemplate;

import com.mongodb.DBCursor;
import com.mongodb.DBObject;

import fr.cirad.mgdb.model.mongo.maintypes.Individual;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongo.subtypes.SampleId;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mongo.MongoTemplateManager;

/**
 * The Class BEDExportHandler.
 */
public class BEDExportHandler extends AbstractMarkerOrientedExportHandler
{

	/** The Constant LOG. */
	static final Logger LOG = Logger.getLogger(BEDExportHandler.class);

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatName()
	 */
	@Override
	public String getExportFormatName() {
		return "BED";
	}

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
	 */
	@Override
	public String getExportFormatDescription() {
		return "Exports data in BED Format. See <a target='_blank' href='http://genome.ucsc.edu/FAQ/FAQformat.html#format1'>http://genome.ucsc.edu/FAQ/FAQformat.html#format1</a> for more details";
	}

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
	 */
	@Override
	public List<String> getStepList() {
		return Arrays.asList(new String[] {"Exporting data to BED format"});
	}

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.markeroriented.AbstractMarkerOrientedExportHandler#exportData(java.io.OutputStream, java.lang.String, java.util.List, fr.cirad.tools.ProgressIndicator, com.mongodb.DBCursor, java.util.Map, int, int, java.util.Map)
	 */
	@Override
	public void exportData(OutputStream outputStream, String sModule, List<SampleId> sampleIDs, ProgressIndicator progress, DBCursor markerCursor, Map<Comparable, Comparable> markerSynonyms, int nMinimumGenotypeQuality, int nMinimumReadDepth, Map<String, InputStream> readyToExportFiles) throws Exception {
		MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
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
			}

		int markerCount = markerCursor.count();

		List<String> selectedIndividualList = new ArrayList<String>();
		for (Individual ind : getIndividualsFromSamples(sModule, sampleIDs))
			selectedIndividualList.add(ind.getId());

		String exportName = sModule + "_" + markerCount + "variants_" + selectedIndividualList.size() + "individuals";
		zos.putNextEntry(new ZipEntry(exportName + ".bed"));

		short nProgress = 0, nPreviousProgress = 0;
		int nChunkSize = Math.min(2000, markerCount), nLoadedMarkerCount = 0;
		while (markerCursor.hasNext())
		{
			int nLoadedMarkerCountInLoop = 0;
			Map<Comparable, String> markerChromosomalPositions = new LinkedHashMap<Comparable, String>();
			boolean fStartingNewChunk = true;
			markerCursor.batchSize(nChunkSize);
			while (markerCursor.hasNext() && (fStartingNewChunk || nLoadedMarkerCountInLoop%nChunkSize != 0)) {
				DBObject exportVariant = markerCursor.next();
				DBObject refPos = (DBObject) exportVariant.get(VariantData.FIELDNAME_REFERENCE_POSITION);
				markerChromosomalPositions.put((Comparable) exportVariant.get("_id"), refPos.get(ReferencePosition.FIELDNAME_SEQUENCE) + ":" + refPos.get(ReferencePosition.FIELDNAME_START_SITE));
				nLoadedMarkerCountInLoop++;
				fStartingNewChunk = false;
			}

			for (Comparable variantId : markerChromosomalPositions.keySet()) // read data and write results into temporary files (one per sample)
			{
				String[] chromAndPos = markerChromosomalPositions.get(variantId).split(":");
				zos.write((chromAndPos[0] + "\t" + (Long.parseLong(chromAndPos[1])-1)  + "\t" + (Long.parseLong(chromAndPos[1])-1) + "\t" + variantId + "\t" + "0" + "\t" + "+").getBytes());
				zos.write((LINE_SEPARATOR).getBytes());
			}

            if (progress.hasAborted())
            	return;

            nLoadedMarkerCount += nLoadedMarkerCountInLoop;
			nProgress = (short) (nLoadedMarkerCount * 100 / markerCount);
			if (nProgress > nPreviousProgress)
			{
				progress.setCurrentStepProgress(nProgress);
				nPreviousProgress = nProgress;
			}
        }

        zos.close();
		progress.setCurrentStepProgress((short) 100);
	}
}