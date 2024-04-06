package at.tugraz.genome.vos;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import at.tugraz.genome.maspectras.utils.Calculator;
import at.tugraz.genome.parsers.MSDialResultsCombiner;

/**
 * @author Juergen Hartler
 */
public class MSDialCombinedEntry
{
	
	//private int totalEntries_;
	private List<String> fileIds_;
	private Hashtable<String,Vector<MSDialEntry>> entries_;
	private Hashtable<String,MSDialEntry> highestScoring_;
	
  private String dialClassName_;
  private String dialMs1Name_;
  private String dialMs2Name_;
  private String ldaClassName_;
  private String ldaMs1Name_;
  private String ldaMs2Name_;
  private String alexClassName_;
  private String alexMs1Name_;
  private String alexMs2Name_;
  private String detectedRts_;
  private String groupingRt_;
  private String adduct_;

	  private double mz_;
  private float totalScoreAvg_;
  private float totalScoreMax_;
  private String detectedTotalScoresMax_;
  private float dotProductAvg_;
  private String detectedDotProductMax_;
  private float weightedDotProductAvg_;
  private String detectedWeightedDotProductMax_;
  private float reverseDotProductAvg_;
  private String detectedreverseDotProductMax_;
  
  private String tp_;
  private String comment_;
  
  private static HashMap<String, String> alexAdductLookup_ = new HashMap<String, String>(){
    private static final long serialVersionUID = -876010361435525395L;

    {
      put("H", "+H+");
      put("Na", "+Na+");
      put("NH4", "+NH4+");
      put("-H", "-H+");
      put("HCOO", "+HCOO-");
    }
  };
  
	
	public MSDialCombinedEntry(Vector<MSDialEntry> entries) {
		this.addEntriesToHash(entries);
		//this.entries_ = entries;
		this.setCommonSettings();
		this.calculateCombinedValues();
  this.tp_ = "false";
  this.comment_ = null;
	}
	
	private void addEntriesToHash(Vector<MSDialEntry> entries) {
		//this.totalEntries_ = entries.size();
		this.entries_ = new Hashtable<String,Vector<MSDialEntry>>();
		this.highestScoring_ = new Hashtable<String,MSDialEntry>();
		fileIds_ = new ArrayList<String>();
		Set<String> ids = new HashSet<String>();
		for (MSDialEntry entry : entries) {
			if (!ids.contains(entry.getFileId())) {
				ids.add(entry.getFileId());
				fileIds_.add(entry.getFileId());
				this.entries_.put(entry.getFileId(), new Vector<MSDialEntry>());
			}
			this.entries_.get(entry.getFileId()).add(entry);
		}
		Collections.sort(fileIds_);
		for (String fileId : this.fileIds_) {
			MSDialEntry highest = null;
			float highestScore = Float.NEGATIVE_INFINITY;
			for (MSDialEntry entry : this.entries_.get(fileId)) {
				if (entry.getScore()>highestScore) {
					highest = entry;
					highestScore = entry.getScore();
				}
			}
			highestScoring_.put(fileId, highest);
		}
		
	}
	
	private void setCommonSettings() {
		MSDialEntry entry = entries_.values().iterator().next().iterator().next();
		this.dialClassName_ = entry.getDialClassName();
		this.dialMs1Name_ = entry.getDialMs1Name();
		this.dialMs2Name_ = entry.getDialMs2Name();
		this.ldaClassName_ = entry.getLdaClassName();
		this.ldaMs1Name_ = entry.getLdaMs1Name();
		this.ldaMs2Name_ = entry.getLdaMs2Name();
		this.alexClassName_ = entry.getAlexClassName();
		this.alexMs1Name_ = entry.getAlexMs1Name();
		this.alexMs2Name_ = entry.getAlexMs2Name();
		this.groupingRt_ = entry.getGroupingRt();
		this.adduct_ = entry.getAdduct();
	}
	
	private void calculateCombinedValues() {
		totalScoreMax_ = Float.NEGATIVE_INFINITY;
		double[] mzs = new double[this.highestScoring_.size()];
		float[] totalScores = new float[this.highestScoring_.size()];
		float[] dotProducts = new float[this.highestScoring_.size()];
		float[] weightedDotProducts = new float[this.highestScoring_.size()];
		float[] reverseDotProducts = new float[this.highestScoring_.size()];
		StringBuilder rts = new StringBuilder();
		StringBuilder totalSs = new StringBuilder();
		StringBuilder dotPs = new StringBuilder();
		StringBuilder weightedDotPs = new StringBuilder();
		StringBuilder reverseDotPs = new StringBuilder();
		for (int i=0; i!=fileIds_.size(); i++) {
			String fileId = fileIds_.get(i);
			if (rts.length()>0) {
				rts.append("|");
				totalSs.append("|");
				dotPs.append("|");
				weightedDotPs.append("|");
				reverseDotPs.append("|");
			}
			rts.append(fileId+"=");
			totalSs.append(fileId+"=");
			dotPs.append(fileId+"=");
			weightedDotPs.append(fileId+"=");
			reverseDotPs.append(fileId+"=");
			int count = 0;
			for (MSDialEntry entry : entries_.get(fileId)) {
				if (count>0) {
					rts.append(";");
					totalSs.append(";");
					dotPs.append(";");
					weightedDotPs.append(";");
					reverseDotPs.append(";");
				}
				rts.append(MSDialResultsCombiner.createUniqueClusterRtString(entry.getRt()));
				totalSs.append(Calculator.roundFloat(entry.getScore(),2));
				dotPs.append(Calculator.roundFloat(entry.getDotProduct(),2));
				weightedDotPs.append(Calculator.roundFloat(entry.getWeightedDotProduct(),2));
				reverseDotPs.append(Calculator.roundFloat(entry.getReverseDotProduct(),2));
				count++;
			}
			MSDialEntry highest = this.highestScoring_.get(fileId);
			mzs[i] = highest.getMz();
			if (highest.getScore()>totalScoreMax_)
				totalScoreMax_ = highest.getScore();
			totalScores[i] = highest.getScore();
			dotProducts[i] = highest.getDotProduct();
			weightedDotProducts[i] = highest.getWeightedDotProduct();
			reverseDotProducts[i] = highest.getReverseDotProduct();
		}
		this.detectedRts_ = rts.toString();
this.detectedTotalScoresMax_ = totalSs.toString();
this.detectedDotProductMax_ = dotPs.toString();
this.detectedWeightedDotProductMax_ = weightedDotPs.toString();
this.detectedreverseDotProductMax_ = reverseDotPs.toString();
		this.mz_ = Calculator.mean(mzs);
		this.totalScoreAvg_ = Calculator.mean(totalScores);
		this.dotProductAvg_ = Calculator.mean(dotProducts);
		this.weightedDotProductAvg_ = Calculator.mean(weightedDotProducts);
		this.reverseDotProductAvg_ = Calculator.mean(reverseDotProducts);
		
//		if (this.detectedRts_.contains(";"))
//			System.out.println(alexMs2Name_+": "+this.detectedRts_);
	}

	public String getDialClassName()
	{
		return dialClassName_;
	}

	public String getDialMs1Name()
	{
		return dialMs1Name_;
	}

	public String getDialMs2Name()
	{
		return dialMs2Name_;
	}

	public String getLdaClassName()
	{
		return ldaClassName_;
	}

	public String getLdaMs1Name()
	{
		return ldaMs1Name_;
	}

	public String getLdaMs2Name()
	{
		return ldaMs2Name_;
	}

	public String getAlexClassName()
	{
		return alexClassName_;
	}

	public String getAlexMs1Name()
	{
		return alexMs1Name_;
	}

	public String getAlexMs2Name()
	{
		return alexMs2Name_;
	}

	public String getDetectedRts()
	{
		return detectedRts_;
	}

	public String getGroupingRt()
	{
		return groupingRt_;
	}

	public String getAdduct()
	{
		return adduct_;
	}

  public String getAlexAdduct(){
    if (alexAdductLookup_.containsKey(adduct_))
      return alexAdductLookup_.get(adduct_);
    else
      return adduct_;
  }


	public double getMz()
	{
		return mz_;
	}

	public float getTotalScoreAvg()
	{
		return totalScoreAvg_;
	}

	public float getTotalScoreMax()
	{
		return totalScoreMax_;
	}

	public String getDetectedTotalScoresMax()
	{
		return detectedTotalScoresMax_;
	}

	public float getDotProductAvg()
	{
		return dotProductAvg_;
	}

	public String getDetectedDotProductMax()
	{
		return detectedDotProductMax_;
	}

	public float getWeightedDotProductAvg()
	{
		return weightedDotProductAvg_;
	}

	public String getDetectedWeightedDotProductMax()
	{
		return detectedWeightedDotProductMax_;
	}

	public float getReverseDotProductAvg()
	{
		return reverseDotProductAvg_;
	}

	public String getDetectedReverseDotProductMax()
	{
		return detectedreverseDotProductMax_;
	}

  public String getTp()
  {
    return tp_;
  }

  public void setTp(String tp)
  {
    this.tp_ = tp;
  }

  public String getComment()
  {
    return comment_;
  }

  public void setComment_(String comment)
  {
    this.comment_ = comment;
  }
	
  
}
