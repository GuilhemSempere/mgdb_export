package fr.cirad.mgdb.exporting.tools;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.log4j.Logger;
import org.springframework.data.mongodb.core.MongoTemplate;

import com.mongodb.DBCursor;
import com.mongodb.DBObject;

import fr.cirad.mgdb.exporting.markeroriented.AbstractMarkerOrientedExportHandler;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
import fr.cirad.mgdb.model.mongo.subtypes.SampleId;
import fr.cirad.mgdb.model.mongodao.MgdbDao;
import fr.cirad.tools.ProgressIndicator;

public class AsyncExportTool {
	
	private static final Logger LOG = Logger.getLogger(AbstractMarkerOrientedExportHandler.class);
	
	/** The Constant NUMBER_OF_SIMULTANEOUS_QUERY_THREADS. */
	static final int NUMBER_OF_SIMULTANEOUS_QUERY_THREADS = 5;
	
	static final int CHUNK_STATUS_BLANK = 0;
	static final int CHUNK_STATUS_QUERYING = 1;
	static final int CHUNK_STATUS_WAITING = 2;
	static final int CHUNK_STATUS_PROCESSED = 3;

	private DBCursor markerCursor;
	private MongoTemplate mongoTemplate;
	private int nQueryChunkSize;
	private List<SampleId> sampleIDs;
	private ProgressIndicator progress;
	private AbstractDataOutputHandler<Integer, LinkedHashMap<VariantData, Collection<VariantRunData>>[]> dataOutputHandler;

	private int queryChunkIndex = 0, lastWrittenChunkIndex = -1;
	private AtomicInteger runningThreadCount = new AtomicInteger(0);
	private LinkedHashMap<VariantData, Collection<VariantRunData>>[] variantDataArray;
	private int[] chunkStatusArray;
	private boolean fHasBeenLaunched = false;

	public AsyncExportTool(DBCursor markerCursor, int nTotalMarkerCount, int nQueryChunkSize, MongoTemplate mongoTemplate, List<SampleId> sampleIDs, AbstractDataOutputHandler dataOutputHandler, ProgressIndicator progress)
	{
		this.mongoTemplate = mongoTemplate;
		this.markerCursor = markerCursor;
		this.nQueryChunkSize = nQueryChunkSize;
		this.sampleIDs = sampleIDs;
		this.progress = progress;
		int nNumberOfChunks = (int) Math.ceil((float) nTotalMarkerCount / nQueryChunkSize);
		this.variantDataArray = new LinkedHashMap[nNumberOfChunks];
		this.chunkStatusArray = new int[nNumberOfChunks];
		this.dataOutputHandler = dataOutputHandler;
	}


	public void launch() throws Exception
	{
		if (fHasBeenLaunched)
			throw new Exception("This AsyncExportTool has already been launched!");
		
		
		new Thread()
		{
			public void run()
			{
				getNextChunk();
			}
		}.start();
		
//		while ((progress == null || !progress.hasAborted()) && (chunkStatusArray.length > lastWrittenChunkIndex + 1))
			new Thread()
			{
				public void run()
				{
					tryAndWriteToOutput();
				}
			}.start();
				
	}
	
	protected boolean getNextChunk()
	{
		if (progress != null && progress.hasAborted())
			return false;

		if (!markerCursor.hasNext())
			return false;

		fHasBeenLaunched = true;
		
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

		final int nFinalChunkIndex = queryChunkIndex++;
		Thread t = new Thread()
		{
			public void run()
			{
				try
				{
					runningThreadCount.addAndGet(1);
					chunkStatusArray[nFinalChunkIndex] = CHUNK_STATUS_QUERYING;
					variantDataArray[nFinalChunkIndex] = MgdbDao.getSampleGenotypes(mongoTemplate, sampleIDs, currentMarkers, true, null /*new Sort(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ChromosomalPosition.FIELDNAME_SEQUENCE).and(new Sort(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ChromosomalPosition.FIELDNAME_START_SITE))*/);
//					System.out.println(nFinalChunkIndex+1 + " / " + variantDataArray.length + " -> read : " + Thread.currentThread().getName());
					chunkStatusArray[nFinalChunkIndex] = CHUNK_STATUS_WAITING;
					runningThreadCount.addAndGet(-1);
					if (runningThreadCount.get() <= NUMBER_OF_SIMULTANEOUS_QUERY_THREADS)
						getNextChunk();
				}
				catch (Exception e)
				{
					LOG.error("Error reading asynchronously genotypes for export", e);
					progress.setError("Error reading asynchronously genotypes for export: " + e);
				}
			}
		};
		
		t.start();
		return true;
	}
	
	static public abstract class AbstractDataOutputHandler<T, V> implements Callable<Void> {
		protected T chunkIndex;
	    protected V variantDataArray;

	    public void setChunkIndex(T chunkIndex) {
	        this.chunkIndex = chunkIndex;
	    }
	    
	    public void setvariantDataArray(V variantDataArray) {
	        this.variantDataArray = variantDataArray;
	    }

	    public abstract Void call ();
	}
	
	protected void tryAndWriteToOutput()
	{
//		if (progress != null && progress.hasAborted() || (chunkStatusArray.length <= lastWrittenChunkIndex + 1))
//			return; // stop here

		try
		{
			while ((progress == null || !progress.hasAborted()) && (chunkStatusArray.length > lastWrittenChunkIndex + 1))
			{
				while (chunkStatusArray[lastWrittenChunkIndex + 1] != CHUNK_STATUS_WAITING)
					Thread.sleep(500);
				
				dataOutputHandler.setChunkIndex(lastWrittenChunkIndex + 1);
				dataOutputHandler.setvariantDataArray(variantDataArray);
				dataOutputHandler.call();	// invoke data writing
				chunkStatusArray[++lastWrittenChunkIndex] = CHUNK_STATUS_PROCESSED;
//				System.out.println(lastWrittenChunkIndex+1 + " / " + variantDataArray.length + " -> written : " + Thread.currentThread().getName());
				progress.setCurrentStepProgress((lastWrittenChunkIndex + 1) * 100 / chunkStatusArray.length);
			}
		}
		catch (Throwable t)
		{
			LOG.error("Error writing export output", t);
			progress.setError("Error writing export output: " + t.getMessage());
		}
	}
}
