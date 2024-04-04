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
	
	private HashMap<String,Set<String>> irrelevantClasses_;
	
	//first key: class; second key: MS1 species; set entries: mod
	private HashMap<String,HashMap<String,Set<String>>> analytesForRtGroupingPos_;
	//first key: class; second key: MS1 species; set entries: mod
	private HashMap<String,HashMap<String,Set<String>>> analytesForRtGroupingNeg_;
	//private LinkedHashMap<String,Vector<String>> analyteSequence_;
	
	//private Hashtable<String,Vector<String>> notCoveredAnalytes_; 
	
	public MSDialResultsCombiner(String positiveIonModeResultsDir, String negativeIonModeResultsDir, HashMap<String,Set<String>> relevantClasses/*, LinkedHashMap<String,Vector<String>> analyteSequence*/) {
		this.posDir_ = positiveIonModeResultsDir;
		this.negDir_ = negativeIonModeResultsDir;
		this.relevantClasses_ = relevantClasses;
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

	
	private void combineResults() {
		System.out.println("Now I have to combine the results");
		this.identifyRTGroupingTimes(this.analytesForRtGroupingPos_,this.posResults_);
		this.identifyRTGroupingTimes(this.analytesForRtGroupingNeg_, this.negResults_);
	}
	
	private void identifyRTGroupingTimes(HashMap<String,HashMap<String,Set<String>>> analytesForRtGrouping,
			LinkedHashMap<String,Hashtable<String,Hashtable<String,Hashtable<String,Vector<MSDialEntry>>>>> results) {
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
	private void groupByRetentionTime(List<MSDialEntry> hitsForAnalyteModAllSearches) {
		List<MSDialEntry> list = new ArrayList<MSDialEntry>(hitsForAnalyteModAllSearches);
		Collections.sort(list,new GeneralComparator("at.tugraz.genome.vos.MSDialEntry", "getArea", "java.lang.Float"));
		int clusterId = 0;
		for (int i=(list.size()-1); i!=-1; i--) {
			//here I stopped
		}
	}
	
}
