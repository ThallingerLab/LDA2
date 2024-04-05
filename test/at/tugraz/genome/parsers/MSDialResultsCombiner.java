package at.tugraz.genome.parsers;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import at.tugraz.genome.exception.MSDialException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.maspectras.utils.Calculator;
import at.tugraz.genome.vos.MSDialEntry;
import at.tugraz.genome.voutils.GeneralComparator;

/**
 * @author Juergen Hartler
 */
public class MSDialResultsCombiner
{
	
	private final static String msDialVersion_ = MSDialEntry.MSDIAL_VERSION_4_9;
	
	private String posDir_;
	private String negDir_;
	private LinkedHashMap<String,Hashtable<String,Hashtable<String,Hashtable<String,Vector<MSDialEntry>>>>> posResults_;
	private LinkedHashMap<String,Hashtable<String,Hashtable<String,Hashtable<String,Vector<MSDialEntry>>>>> negResults_;
	
	private HashMap<String,Set<String>> relevantClasses_;
	private double rtGroupingTolerance_;
	
	private HashMap<String,Set<String>> irrelevantClasses_;
	
	//first key: class; second key: MS1 species; set entries: mod
	private HashMap<String,HashMap<String,Set<String>>> analytesForRtGroupingPos_;
	//first key: class; second key: MS1 species; set entries: mod
	private HashMap<String,HashMap<String,Set<String>>> analytesForRtGroupingNeg_;
	//private LinkedHashMap<String,Vector<String>> analyteSequence_;
	
	//private Hashtable<String,Vector<String>> notCoveredAnalytes_; 
	
	public MSDialResultsCombiner(String positiveIonModeResultsDir, String negativeIonModeResultsDir, HashMap<String,Set<String>> relevantClasses, double rtGroupingTolerance/*, LinkedHashMap<String,Vector<String>> analyteSequence*/) {
		this.posDir_ = positiveIonModeResultsDir;
		this.negDir_ = negativeIonModeResultsDir;
		this.relevantClasses_ = relevantClasses;
		this.rtGroupingTolerance_ = rtGroupingTolerance;
		//this.analyteSequence_ = analyteSequence;
	}
	
	public void parseAndCombine() throws MSDialException, LipidCombinameEncodingException {
		this.parseMSDIALResults();
		this.checkMissingMSDIALClasses();
		this.combineResults();
	}
	
	private void parseMSDIALResults() throws MSDialException, LipidCombinameEncodingException {
		posResults_ = new LinkedHashMap<String,Hashtable<String,Hashtable<String,Hashtable<String,Vector<MSDialEntry>>>>>();
		negResults_ = new LinkedHashMap<String,Hashtable<String,Hashtable<String,Hashtable<String,Vector<MSDialEntry>>>>>();
		irrelevantClasses_= new HashMap<String,Set<String>>();
		irrelevantClasses_.put("+", new HashSet<String>());
		irrelevantClasses_.put("-", new HashSet<String>());
		analytesForRtGroupingPos_ = new HashMap<String,HashMap<String,Set<String>>>();
		analytesForRtGroupingNeg_ = new HashMap<String,HashMap<String,Set<String>>>();
		
		File posDir = new File(posDir_);
		File negDir = new File(negDir_);
		if (!posDir.isDirectory())
			throw new MSDialException("MSDialResultsCombiner: You must specify a directory for the positive ion mode data! This is a file: "+posDir_);
		if (!negDir.isDirectory())
			throw new MSDialException("MSDialResultsCombiner: You must specify a directory for the negative ion mode data! This is a file: "+negDir_);
		for (File aFile : posDir.listFiles()) {
			if (aFile.isDirectory() || !aFile.getName().endsWith(".txt"))
				continue;
			String shortName = aFile.getName().substring(0,aFile.getName().indexOf("_"));
			MSDialTxtParser parser = new MSDialTxtParser(aFile.getAbsolutePath(),shortName);
			parser.parse(msDialVersion_);
			Hashtable<String,Hashtable<String,Hashtable<String,Vector<MSDialEntry>>>> results = parser.getStructuredResults();
			Hashtable<String,Hashtable<String,Hashtable<String,Vector<MSDialEntry>>>> relevantResults = new Hashtable<String,Hashtable<String,Hashtable<String,Vector<MSDialEntry>>>>(); 
			posResults_.put(shortName, relevantResults);
			for (String lClass : results.keySet()) {
				if (relevantClasses_.containsKey(lClass) && relevantClasses_.get(lClass).contains("+")) {
					relevantResults.put(lClass, results.get(lClass));
					if (!analytesForRtGroupingPos_.containsKey(lClass))
						analytesForRtGroupingPos_.put(lClass, new HashMap<String,Set<String>>());
					for (String analyte : results.get(lClass).keySet()) {
						if (!analytesForRtGroupingPos_.get(lClass).containsKey(analyte))
							analytesForRtGroupingPos_.get(lClass).put(analyte, new HashSet<String>());
						for (String mod : results.get(lClass).get(analyte).keySet())
							analytesForRtGroupingPos_.get(lClass).get(analyte).add(mod);
					}
				} else {
					irrelevantClasses_.get("+").add(lClass);
				}
			}
		}
		for (File aFile : negDir.listFiles()) {
			if (aFile.isDirectory() || !aFile.getName().endsWith(".txt"))
				continue;
			String shortName = aFile.getName().substring(0,aFile.getName().indexOf("_"));
			MSDialTxtParser parser = new MSDialTxtParser(aFile.getAbsolutePath(),shortName);
			parser.parse(msDialVersion_);
			Hashtable<String,Hashtable<String,Hashtable<String,Vector<MSDialEntry>>>> results = parser.getStructuredResults();
			Hashtable<String,Hashtable<String,Hashtable<String,Vector<MSDialEntry>>>> relevantResults = new Hashtable<String,Hashtable<String,Hashtable<String,Vector<MSDialEntry>>>>(); 
			negResults_.put(shortName, relevantResults);
			for (String lClass : results.keySet()) {
				if (relevantClasses_.containsKey(lClass) && relevantClasses_.get(lClass).contains("-")) {
					relevantResults.put(lClass, results.get(lClass));
					if (!analytesForRtGroupingNeg_.containsKey(lClass))
						analytesForRtGroupingNeg_.put(lClass, new HashMap<String,Set<String>>());
					for (String analyte : results.get(lClass).keySet()) {
						if (!analytesForRtGroupingNeg_.get(lClass).containsKey(analyte))
							analytesForRtGroupingNeg_.get(lClass).put(analyte, new HashSet<String>());
						for (String mod : results.get(lClass).get(analyte).keySet())
							analytesForRtGroupingNeg_.get(lClass).get(analyte).add(mod);
					}

				} else {
					irrelevantClasses_.get("-").add(lClass);
				}
			}
		}		
	}
	
	private void checkMissingMSDIALClasses() {
		
		for (String lClass : irrelevantClasses_.get("+")) {
				System.out.println("The following MSDial class is neglected in + ion mode: "+lClass);
		}
		for (String lClass : irrelevantClasses_.get("-")) {
			System.out.println("The following MSDial class is neglected in - ion mode: "+lClass);
	}

	}

	
	private void combineResults() throws MSDialException {
		System.out.println("Now I have to combine the results");
		this.identifyRTGroupingTimes(this.analytesForRtGroupingPos_,this.posResults_);
		this.identifyRTGroupingTimes(this.analytesForRtGroupingNeg_, this.negResults_);
	}
	
	private void identifyRTGroupingTimes(HashMap<String,HashMap<String,Set<String>>> analytesForRtGrouping,
			LinkedHashMap<String,Hashtable<String,Hashtable<String,Hashtable<String,Vector<MSDialEntry>>>>> results) throws MSDialException {
		for (String lClass : analytesForRtGrouping.keySet()) {
			System.out.println("LipidClass: "+lClass);
			for (String analyte : analytesForRtGrouping.get(lClass).keySet()) {
				for (String mod : analytesForRtGrouping.get(lClass).get(analyte)) {
					List<MSDialEntry> hitsForAnalyteModAllSearches = new ArrayList<MSDialEntry>();
					for (String msRunKey : results.keySet()) {
						Hashtable<String,Hashtable<String,Hashtable<String,Vector<MSDialEntry>>>> oneSearch = results.get(msRunKey);
						if (!oneSearch.containsKey(lClass))
							continue;
						if (!oneSearch.get(lClass).containsKey(analyte))
							continue;
						if (!oneSearch.get(lClass).get(analyte).containsKey(mod))
							continue;
						Vector<MSDialEntry> entries = oneSearch.get(lClass).get(analyte).get(mod);
						if (entries.size()>0)
							hitsForAnalyteModAllSearches.addAll(entries);
					}
					this.groupByRetentionTime(hitsForAnalyteModAllSearches);
					//System.out.println(mod+" - "+analyte+" - "+hitsForAnalyteModAllSearches.size());
				}
			}
		}
	}
	
	@SuppressWarnings("unchecked")
	private void groupByRetentionTime(List<MSDialEntry> hitsForAnalyteModAllSearches) throws MSDialException {
		List<MSDialEntry> list = new ArrayList<MSDialEntry>(hitsForAnalyteModAllSearches);
		Collections.sort(list,new GeneralComparator("at.tugraz.genome.vos.MSDialEntry", "getArea", "java.lang.Float"));
		
    Hashtable<Integer,Hashtable<String,Vector<MSDialEntry>>> valuesInClusters = new  Hashtable<Integer,Hashtable<String,Vector<MSDialEntry>>>();
    //Hashtable<String,Set<String>> usedRts = new Hashtable<String,Set<String>>();
    Hashtable<Integer,Double> rtClusters = new Hashtable<Integer,Double>();
    LinkedHashMap<MSDialEntry,Boolean> hitsInsideAClusterRange = new LinkedHashMap<MSDialEntry,Boolean>();
		for (int i=(list.size()-1); i!=-1; i--) {
			hitsInsideAClusterRange.put(list.get(i), new Boolean(false));
		}
		int clusterId = 0;
    while (areThereHitsOutsideClusterRange(hitsInsideAClusterRange, rtClusters/*, usedRts*/)){
			MSDialEntry entry =  findStrongestHitOutsideExistingClusters(hitsInsideAClusterRange);
			String rt = createUniqueClusterRtString(entry.getRt());
			addClosestPeaksToCluster(clusterId, entry, rt, /*usedRts,*/ valuesInClusters, hitsInsideAClusterRange);
			addAreaWeightedMeanRt(rtClusters,clusterId,valuesInClusters.get(clusterId));
//			if (clusterId>3)
//				System.out.println("ClusterId: "+clusterId+" ; "+entry.getAlexClassName()+" ; "+entry.getAlexMs1Name()+" ; "+entry.getAlexMs2Name()+" ; "+entry.getRt());
			clusterId++;
		}
    addRemainingPeaksToClosestClusters(rtClusters,valuesInClusters,/*usedRts,*/hitsInsideAClusterRange);
    rtClusters = new Hashtable<Integer,Double>();
    for (int cluster : valuesInClusters.keySet()){
      addAreaWeightedMeanRt(rtClusters,cluster,valuesInClusters.get(cluster));
    }
    while (clusterOverlap(rtClusters)){
      int[] overlapIds = detectClosestOverlap(rtClusters);
      uniteTwoClusters(overlapIds[0],overlapIds[1],valuesInClusters);
      rtClusters = new Hashtable<Integer,Double>();
      for (int cluster : valuesInClusters.keySet()){
        addAreaWeightedMeanRt(rtClusters,cluster,valuesInClusters.get(cluster));
      }
    }
    if (clusterOverlap(rtClusters)) {
    	System.out.println("I do have overlapping clusters");
    }
    //grouping is done; now generate unique RT strings
    Hashtable<Integer,String> rtClusterStrings = new Hashtable<Integer,String>();
    for (Integer cId : rtClusters.keySet()) {
    	rtClusterStrings.put(cId, createUniqueClusterRtString((rtClusters.get(cId)).floatValue()));
    }
    //now add the grouping RT to every MS-Dial Entry
    for (int i=(list.size()-1); i!=-1; i--) {
    	MSDialEntry entry =  list.get(i);
    	entry.setGroupingRt(rtClusterStrings.get(entry.getRtClusterId()));
			System.out.println("ClusterId: "+entry.getRtClusterId()+" ; "+entry.getFileId()+" ; "+entry.getAlexClassName()+" ; "+entry.getAlexMs1Name()+" ; "+entry.getAlexMs2Name()+" ; "+entry.getRt()+": "+entry.getGroupingRt());

    }

	}
	
  private boolean areThereHitsOutsideClusterRange(LinkedHashMap<MSDialEntry,Boolean> hits, Hashtable<Integer,Double> rtClusters/*, Hashtable<String,Set<String>> usedRts*/){
    boolean anyFound = false;
    boolean outsideDetected = false;
    for (MSDialEntry entry : hits.keySet()) {
      anyFound = true;
      if (hits.get(entry))
      	continue;
      double rtHit = (double)entry.getRt();
      boolean insideACluster = false;
      for (double rt : rtClusters.values()){
      	if (isWithinRtGroupingBoundaries(rtHit,rt)){
      		insideACluster = true;
      		hits.put(entry, true);
      		break;
      	}
      }
      if (!insideACluster){
      	outsideDetected = true;
      }
      if (outsideDetected) break;
    }
    if (!anyFound)
      return false;
    else
      return outsideDetected;
  }
    
  private void addClosestPeaksToCluster(int clusterId, MSDialEntry strongestHit,
      String strongestRt, /*Vector<String> fileNames, Hashtable<String,Set<String>> usedRts,*/
      Hashtable<Integer,Hashtable<String,Vector<MSDialEntry>>> clusters,
      LinkedHashMap<MSDialEntry,Boolean> hits){
    //ResultAreaVO strongestHit = resultsInHash.get(strongestFile).get(groupName).get(molName).get(strongestRt);
    double strongRt = Double.parseDouble(strongestRt);
    
    //first key: fileName; second key rt
    Hashtable<String,Hashtable<String,MSDialEntry>> closest = new Hashtable<String,Hashtable<String,MSDialEntry>>();
    //find hits in other searches which are within the tolerance and closest to strongest one
    for (MSDialEntry entry : hits.keySet()) {
    	if (!isWithinRtGroupingBoundaries(entry.getRt(),strongRt))
    		continue;
    	hits.put(entry, true);
    	if (entry.getFileId().equals(strongestHit.getFileId()) || entry.getRt()==strongRt)
    		continue;
    	if (closest.containsKey(entry.getFileId()) && Math.abs(entry.getRt()-strongRt)<Math.abs(Double.valueOf(closest.get(entry.getFileId()).keySet().iterator().next())-strongRt))
    		continue;
    	Hashtable<String,MSDialEntry> closestHash = new Hashtable<String,MSDialEntry>();
    	String closestRt = createUniqueClusterRtString(entry.getRt());
      closestHash.put(closestRt, entry);
      closest.put(entry.getFileId(), closestHash);
    }
    Hashtable<String,Vector<MSDialEntry>> cluster = new Hashtable<String,Vector<MSDialEntry>>();
    ////Set<String> rts = new HashSet<String>();
    Vector<MSDialEntry> oneFile = new Vector<MSDialEntry>();
    oneFile.add(strongestHit);
    ////rts.add(strongestRt);
    cluster.put(strongestHit.getFileId(), oneFile);
    ////usedRts.put(strongestHit.getFileId(), rts);
    //TODO: I am not sure wheter I should do that already at this place
    strongestHit.setRtClusterId(clusterId);
    //strongestHit.setGroupingRt(strongestRt);
    //calculate the median of the closest peaks
    if (closest.size()>0){
      float[] differenceToStrongest = new float[closest.size()];
      int count = 0;
      for (String fileName : closest.keySet()){
        differenceToStrongest[count] = Math.abs(Float.parseFloat(closest.get(fileName).keySet().iterator().next())-(float)strongRt);
        count++;
      }
      float medianDifference = Calculator.median(differenceToStrongest);
      //add all hits to the cluster which are within the range of two times the median difference
      for (String fileName : closest.keySet()){
        String rt = closest.get(fileName).keySet().iterator().next();
        float diff = Math.abs(Float.parseFloat(rt)-(float)strongRt);
        if (diff>1.2d*medianDifference)
          continue;
        ////rts = new HashSet<String>();
        oneFile = new Vector<MSDialEntry>();
        oneFile.add(closest.get(fileName).get(rt));
      //TODO: I am not sure wheter I should do that already at this place
        closest.get(fileName).get(rt).setRtClusterId(clusterId);
        //closest.get(fileName).get(rt).setGroupingRt(strongestRt);
        ////rts.add(rt);
        cluster.put(fileName, oneFile);
        ////usedRts.put(fileName, rts);
      }
    }
    clusters.put(clusterId, cluster);
  }
  
  private void addRemainingPeaksToClosestClusters(Hashtable<Integer,Double> rtClusters,
  		Hashtable<Integer,Hashtable<String,Vector<MSDialEntry>>> valuesInClusters, //Hashtable<String,Set<String>> usedRts,
  		LinkedHashMap<MSDialEntry,Boolean> hitsInsideAClusterRange/*,String groupName, String molName, Vector<String> fileNames*/){
    double rtCluster;
    double diff;
    for (MSDialEntry entry : hitsInsideAClusterRange.keySet()) {
    	if (entry.getRtClusterId()!=null)
    		continue;
    	double rtHit = (double)entry.getRt();
    	int bestClusterId = -1;
    	double smallestDiff  = Double.MAX_VALUE;
    	for (Integer clusterId : rtClusters.keySet()){
    		rtCluster = rtClusters.get(clusterId);
    		diff = Math.abs(rtCluster-rtHit);
    		if (diff<smallestDiff){
    			smallestDiff = diff;
    			bestClusterId = clusterId;
    		}
    	}
    	//add MSDialEntry to hash
    	Hashtable<String,Vector<MSDialEntry>> oneCluster = valuesInClusters.get(bestClusterId);
    	Vector<MSDialEntry> inCluster = new Vector<MSDialEntry>();
    	if (oneCluster.containsKey(entry.getFileId())) inCluster = oneCluster.get(entry.getFileId());
    	inCluster.add(entry);
    	entry.setRtClusterId(bestClusterId);
    	oneCluster.put(entry.getFileId(), inCluster);
    	valuesInClusters.put(bestClusterId,oneCluster);
    }
  }
  
  private void addAreaWeightedMeanRt(Hashtable<Integer,Double> rtClusters, int clusterId, Hashtable<String,Vector<MSDialEntry>> cluster){
    double totalArea = 0d;
    double rtTimesArea = 0d;
    double rt;
    double area;
    for (Vector<MSDialEntry> vos : cluster.values()){
      for (MSDialEntry vo : vos){
        rt = (double)vo.getRt();
        area = (double)vo.getArea();
        totalArea += area;
        rtTimesArea += area*rt;
      }
    }
    rtClusters.put(clusterId, rtTimesArea/totalArea);
  }
  
  private boolean clusterOverlap(Hashtable<Integer,Double> rtClusters){
    for (int i=0; i!=rtClusters.size(); i++){
      double firstRt = rtClusters.get(i);
      for (int j=(i+1); j<(rtClusters.size()); j++){
        double secondRt = rtClusters.get(j);
        if (isWithinRtGroupingBoundaries(firstRt, secondRt)){
          return true;
        }
      }
    }
    return false;
  }

  private int[] detectClosestOverlap(Hashtable<Integer,Double> rtClusters){
    double lowestOverlap = this.rtGroupingTolerance_;
    int[] clusterIds = new int[2];
    for (int i=0; i!=rtClusters.size(); i++){
      double firstRt = rtClusters.get(i);
      for (int j=(i+1); j<(rtClusters.size()); j++){
        double secondRt = rtClusters.get(j);
        if (!isWithinRtGroupingBoundaries(firstRt, secondRt))
          continue;
        double overlap = Math.abs(firstRt-secondRt);
        if (overlap>lowestOverlap)
          continue;
        clusterIds[0] = i;
        clusterIds[1] = j;
      }
    }
    return clusterIds;
  }
  
  private void uniteTwoClusters(int id1, int id2, Hashtable<Integer,Hashtable<String,Vector<MSDialEntry>>> valuesInClusters){
    int lower = id1;
    int upper = id2;
    if (lower>upper){
      lower = id2;
      upper = id1;
    }
    Hashtable<String,Vector<MSDialEntry>> strongerCluster = valuesInClusters.get(lower);
    Hashtable<String,Vector<MSDialEntry>> weakerCluster = valuesInClusters.get(upper);
    for (String fileName : weakerCluster.keySet()){
    	Vector<MSDialEntry> toAdd = new Vector<MSDialEntry>();
      if (strongerCluster.containsKey(fileName)) toAdd = strongerCluster.get(fileName);
      Vector<MSDialEntry> toBeAdded = weakerCluster.get(fileName);
      for (MSDialEntry addIt : toBeAdded)
      	addIt.setRtClusterId(lower);
      toAdd.addAll(toBeAdded);
    }
    int count = upper;
    //reorganize the clusterIds
    while ((count+1)<valuesInClusters.size()){
      valuesInClusters.put(count, valuesInClusters.get(count+1));
      for (Vector<MSDialEntry> entries : valuesInClusters.get(count).values()) {
      	for (MSDialEntry entry : entries)
      		entry.setRtClusterId(count);
      }
      count++;
    }
    valuesInClusters.remove(valuesInClusters.size()-1);
  }
  
  private MSDialEntry findStrongestHitOutsideExistingClusters(LinkedHashMap<MSDialEntry,Boolean> hits)throws MSDialException {
  	for (MSDialEntry entry : hits.keySet()) {
  		if (!hits.get(entry))
  			return entry;
  	}
  	throw new MSDialException("There is something wrong with the method MSDialResultsCombiner.areThereHitsOutsideClusterRange. Please contact the developer.");
  }
  
  private boolean isWithinRtGroupingBoundaries(double rt, double refTime){
    return StaticUtils.isWithinTolerance(rtGroupingTolerance_, refTime, rt);
  }
	

  
  private String createUniqueClusterRtString(float rt) {
  	return Calculator.FormatNumberToString(Calculator.roundFloat(rt,3),3d);
  }
}
