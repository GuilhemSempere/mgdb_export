/*******************************************************************************
 * MGDB Export - Mongo Genotype DataBase, export handlers
 * Copyright (C) 2016 <South Green>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License, version 3 as
 * published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * See <http://www.gnu.org/licenses/gpl-3.0.html> for details about
 * GNU General Public License V3.
 *******************************************************************************/
package fr.cirad.mgdb.exporting.markeroriented;

import htsjdk.variant.variantcontext.VariantContext.Type;

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
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.springframework.data.mongodb.core.MongoTemplate;

import com.mongodb.DBCursor;
import com.mongodb.DBObject;

import fr.cirad.mgdb.model.mongo.maintypes.Individual;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongo.subtypes.VariantRunData;
import fr.cirad.mgdb.model.mongo.subtypes.SampleGenotype;
import fr.cirad.mgdb.model.mongo.subtypes.SampleId;
import fr.cirad.mgdb.model.mongodao.MgdbDao;
import fr.cirad.tools.Helper;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mongo.MongoTemplateManager;

/**
 * The Class GFFExportHandler.
 */
public class GFFExportHandler extends AbstractMarkerOrientedExportHandler {

	/** The Constant LOG. */
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
		return "Exports data in GFF3 Format based on Sequence Ontology. See <a target='_blank' href='http://rice.bio.indiana.edu:7082/annot/gff3.html'>http://rice.bio.indiana.edu:7082/annot/gff3.html</a> and <a target='_blank' href='http://www.sequenceontology.org/resources/gff3.html'>http://www.sequenceontology.org/resources/gff3.html</a>";
	}

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.markeroriented.AbstractMarkerOrientedExportHandler#exportData(java.io.OutputStream, java.lang.String, java.util.List, fr.cirad.tools.ProgressIndicator, com.mongodb.DBCursor, java.util.Map, int, int, java.util.Map)
	 */
	@Override
	public void exportData(OutputStream outputStream, String sModule, List<SampleId> sampleIDs, ProgressIndicator progress, DBCursor markerCursor, Map<Comparable, Comparable> markerSynonyms, int nMinimumGenotypeQuality, int nMinimumReadDepth, Map<String, InputStream> readyToExportFiles) throws Exception
	{
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

		File warningFile = File.createTempFile("export_warnings_", "");
		FileWriter warningFileWriter = new FileWriter(warningFile);

		int markerCount = markerCursor.count();

		List<Individual> individuals = getIndividualsFromSamples(sModule, sampleIDs);
		ArrayList<String> individualList = new ArrayList<String>();
		for (int i = 0; i < sampleIDs.size(); i++) {
			Individual individual = individuals.get(i);
			if (!individualList.contains(individual.getId())) {
				individualList.add(individual.getId());
			}
		}

		String exportName = sModule + "_" + markerCount + "variants_" + individualList.size() + "individuals";
		zos.putNextEntry(new ZipEntry(exportName + ".gff3"));
		String header = "##gff-version 3" + LINE_SEPARATOR;
		zos.write(header.getBytes());

		TreeMap<String, String> typeToOntology = new TreeMap<String, String>();
		typeToOntology.put(Type.SNP.toString(), "SO:0000694");
		typeToOntology.put(Type.INDEL.toString(), "SO:1000032");
		typeToOntology.put(Type.MIXED.toString(), "SO:0001059");
		typeToOntology.put(Type.SYMBOLIC.toString(), "SO:0000109");
		typeToOntology.put(Type.MNP.toString(), "SO:0001059");

		int avgObjSize = (Integer) mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantRunData.class)).getStats().get("avgObjSize");
		int nChunkSize = nMaxChunkSizeInMb*1024*1024 / avgObjSize;
		short nProgress = 0, nPreviousProgress = 0;
		long nLoadedMarkerCount = 0;

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

			List<Comparable> currentMarkers = new ArrayList<Comparable>(markerChromosomalPositions.keySet());
			LinkedHashMap<VariantData, Collection<VariantRunData>> variantsAndRuns = MgdbDao.getSampleGenotypes(mongoTemplate, sampleIDs, currentMarkers, true, null /*new Sort(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ChromosomalPosition.FIELDNAME_SEQUENCE).and(new Sort(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ChromosomalPosition.FIELDNAME_START_SITE))*/);	// query mongo db for matching genotypes
			for (VariantData variant : variantsAndRuns.keySet()) // read data and write results into temporary files (one per sample)
			{
				Comparable variantId = variant.getId();
				List<String> variantDataOrigin = new ArrayList<String>();

				Map<String,Integer> gqValueForSampleId = new LinkedHashMap<String, Integer>();
				Map<String,Integer> dpValueForSampleId = new LinkedHashMap<String, Integer>();
				Map<String, List<String>> individualGenotypes = new LinkedHashMap<String, List<String>>();
				List<String> chromAndPos = Helper.split(markerChromosomalPositions.get(variantId), ":");
				if (chromAndPos.size() == 0)
					LOG.warn("Chromosomal position not found for marker " + variantId);
				// LOG.debug(marker + "\t" + (chromAndPos.length == 0 ? "0" : chromAndPos[0]) + "\t" + 0 + "\t" + (chromAndPos.length == 0 ? 0l : Long.parseLong(chromAndPos[1])) + LINE_SEPARATOR);
	        	if (markerSynonyms != null)
	        	{
	        		Comparable syn = markerSynonyms.get(variantId);
	        		if (syn != null)
	        			variantId = syn;
	        	}

				Collection<VariantRunData> runs = variantsAndRuns.get(variant);
				if (runs != null)
					for (VariantRunData run : runs)
						for (Integer sampleIndex : run.getSampleGenotypes().keySet())
						{
							SampleGenotype sampleGenotype = run.getSampleGenotypes().get(sampleIndex);
							String individualId = individuals.get(sampleIDs.indexOf(new SampleId(run.getId().getProjectId(), sampleIndex))).getId();

							Integer gq = null;
							try
							{
								gq = (Integer) sampleGenotype.getAdditionalInfo().get(VariantData.GT_FIELD_GQ);
							}
							catch (Exception ignored)
							{}
							if (gq != null && gq < nMinimumGenotypeQuality)
								continue;

							Integer dp = null;
							try
							{
								dp = (Integer) sampleGenotype.getAdditionalInfo().get(VariantData.GT_FIELD_DP);
							}
							catch (Exception ignored)
							{}
							if (dp != null && dp < nMinimumReadDepth)
								continue;

							String gtCode = sampleGenotype.getCode();
							List<String> storedIndividualGenotypes = individualGenotypes.get(individualId);
							if (storedIndividualGenotypes == null) {
								storedIndividualGenotypes = new ArrayList<String>();
								individualGenotypes.put(individualId, storedIndividualGenotypes);
							}
							storedIndividualGenotypes.add(gtCode);
						}

				zos.write((chromAndPos.get(0) + "\t" + StringUtils.join(variantDataOrigin, ";") /*source*/+ "\t" + typeToOntology.get(variant.getType()) + "\t" + Long.parseLong(chromAndPos.get(1)) + "\t" + Long.parseLong(chromAndPos.get(1)) + "\t" + "." + "\t" + "+" + "\t" + "." + "\t").getBytes());
            	Comparable syn = markerSynonyms == null ? null : markerSynonyms.get(variant.getId());
				zos.write(("ID=" + variant.getId() + ";" + (syn != null ? "Name=" + syn + ";" : "") + "alleles=" + StringUtils.join(variant.getKnownAlleleList(), "/") + ";" + "refallele=" + variant.getKnownAlleleList().get(0) + ";").getBytes());

				for (int j=0; j<individualList.size(); j++ /* we use this list because it has the proper ordering*/)
				{

					NumberFormat nf = NumberFormat.getInstance(Locale.US);
					nf.setMaximumFractionDigits(4);
					HashMap<String, Integer> compt1 = new HashMap<String, Integer>();
					int highestGenotypeCount = 0;
					int sum = 0;

					String individualId = individualList.get(j);
					List<String> genotypes = individualGenotypes.get(individualId);
					HashMap<Object, Integer> genotypeCounts = new HashMap<Object, Integer>(); // will help us to keep track of missing genotypes

					String mostFrequentGenotype = null;
					if (genotypes != null)
						for (String genotype : genotypes) {
							if (genotype.length() == 0)
								continue; /* skip missing genotypes */

							int count = 0;
							for(String t : variant.getAllelesFromGenotypeCode(genotype))
							{
								for (String t1 : variant.getKnownAlleleList()) {
									if (t.equals(t1) && !(compt1.containsKey(t1))){
										count++;
										compt1.put(t1, count);
									}
									else if(t.equals(t1) && compt1.containsKey(t1)){
										if (compt1.get(t1) != 0){
											count++;
											compt1.put(t1, count);
										}else
											compt1.put(t1, count);
									}
									else if (!(compt1.containsKey(t1))){
										compt1.put(t1, 0);
									}
								}
							}
							for(int countValue : compt1.values()){
								sum += countValue;
							}

							int gtCount = 1 + MgdbDao.getCountForKey(genotypeCounts, genotype);
							if (gtCount > highestGenotypeCount) {
								highestGenotypeCount = gtCount;
								mostFrequentGenotype = genotype;
							}
							genotypeCounts.put(genotype, gtCount);
						}

					List<String> alleles = mostFrequentGenotype == null ? new ArrayList<String>() : variant.getAllelesFromGenotypeCode(mostFrequentGenotype);

					if (alleles.size() != 0) {
						zos.write(("acounts=" + individualId+":").getBytes());

						for (String knowAllelesCompt : compt1.keySet()) {
							zos.write((knowAllelesCompt + " " + nf.format(compt1.get(knowAllelesCompt)/(float)sum)+  " " + compt1.get(knowAllelesCompt) + " ").getBytes()) ;
						}
						zos.write((alleles.size() + ";").getBytes());
					}
					if (genotypeCounts.size() > 1)
					{
						Comparable sVariantId = markerSynonyms != null ? markerSynonyms.get(variant.getId()) : variant.getId();
						warningFileWriter.write("- Dissimilar genotypes found for variant " + (sVariantId == null ? variant.getId() : sVariantId) + ", individual " + individualId + ". Exporting most frequent: " + StringUtils.join(alleles, ",") + "\n");
					}
				}
				zos.write((LINE_SEPARATOR).getBytes());
			}

            if (progress.hasAborted())
            	return;

            nLoadedMarkerCount += nLoadedMarkerCountInLoop;
			nProgress = (short) (nLoadedMarkerCount * 100 / markerCount);
			if (nProgress > nPreviousProgress)
			{
//				if (nProgress%5 == 0)
//					LOG.info("========================= exportData: " + nProgress + "% =========================" + (System.currentTimeMillis() - before)/1000 + "s");
				progress.setCurrentStepProgress(nProgress);
				nPreviousProgress = nProgress;
			}
        }

		warningFileWriter.close();
        if (warningFile.length() > 0)
        {
	        zos.putNextEntry(new ZipEntry(exportName + "-REMARKS.txt"));
	        int nWarningCount = 0;
	        BufferedReader in = new BufferedReader(new FileReader(warningFile));
	        String sLine;
	        while ((sLine = in.readLine()) != null)
			{
	        	zos.write((sLine + "\n").getBytes());
	        	in.readLine();
	        	nWarningCount++;
			}
			LOG.info("Number of Warnings for export (" + exportName + "): " + nWarningCount);
			in.close();
        }
        warningFile.delete();

        zos.close();
		progress.setCurrentStepProgress((short) 100);
	}


/* (non-Javadoc)
 * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
 */
@Override
	public List<String> getStepList() {
		return Arrays.asList(new String[] {"Exporting data to GFF3 format"});
	}
}