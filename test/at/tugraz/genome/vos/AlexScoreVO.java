package at.tugraz.genome.vos;

public class AlexScoreVO
{
  
  private String lipidClass_;
  private String lipidSpecies_;
  private String molSpecies_;
  private short identificationType_;
  private double score_;
  private double scoreAllFrags_;
  private double scoreUncorrected_;
  private String adduct_;
  private boolean polarity_;
  private String rtGroup_;
  private double prophetScore_;
  private String tp_;
  private String comment_;
  private int maxFragments_;
  private int maxAllFrags_;
  private String ambiguousPartners_;
  
  public final static short IDENT_TYPE_SPECIES = 0;
  public final static short IDENT_TYPE_MOL_SPECIES = 1;
  
  
  public AlexScoreVO(String lipidClass, String lipidSpecies,
      String molSpecies, short identificationType, double score, double scoreAllFrags,
      String adduct, boolean polarity, String rtGroup,
      String comment, int maxFragments,
      int maxAllFrags, double scoreUncorrected, String ambiguousPartners)
  {
    super();
    this.lipidClass_ = lipidClass;
    this.lipidSpecies_ = lipidSpecies;
    this.molSpecies_ = molSpecies;
    this.identificationType_ = identificationType;
    this.score_ = score;
    this.scoreAllFrags_ = scoreAllFrags;
    this.adduct_ = adduct;
    this.polarity_ = polarity;
    this.rtGroup_ = rtGroup;
    this.prophetScore_ = 0.0d;
    this.tp_ = "";
    this.comment_ = comment;
    this.maxFragments_ = maxFragments;
    this.maxAllFrags_ = maxAllFrags;
    this.scoreUncorrected_ = scoreUncorrected;
    this.ambiguousPartners_ = ambiguousPartners;
  }
  
  public String getLipidClass()
  {
    return lipidClass_;
  }
  public String getLipidSpecies()
  {
    return lipidSpecies_;
  }
  public String getMolSpecies()
  {
    return molSpecies_;
  }
  public short getIdentificationType()
  {
    return identificationType_;
  }
  public double getScore()
  {
    return score_;
  }
  public double getScoreAllFrags()
  {
    return scoreAllFrags_;
  }
  
  public String getComment()
  {
    return comment_;
  }

  public void setComment(String comment)
  {
    this.comment_ = comment;
  }

  public void setScore(double score_)
  {
    this.score_ = score_;
  }
  public void setScoreAllFrags(double scoreAllFrags)
  {
    this.scoreAllFrags_ = scoreAllFrags;
  }
  
  
  public double getScoreUncorrected()
	{
		return scoreUncorrected_;
	}

	public String getAdduct()
  {
    return adduct_;
  }
  public boolean isPositive()
  {
    return polarity_;
  }
  public String getRtGroup()
  {
    return rtGroup_;
  }
  
  public String getId() {
    return this.getMolSpecies()+" "+this.getAdduct();
  }
  
  
  
   public double getProphetScore()
  {
     return prophetScore_;
  }

  public void setProphetScore(double prophetScore)
  {
    this.prophetScore_ = prophetScore;
  }

  public String toString() {
     return "Score: "+this.score_+"; ScoreAllFrags: "+this.scoreAllFrags_+"; ScoreUncorrected: "+this.scoreUncorrected_+"; Class: "+this.lipidClass_+"; Species: "+this.lipidSpecies_+";"+this.molSpecies_+"; Adduct: "+this.adduct_+
         "; Ident-level: "+this.identificationType_+"; Positive-Mode: "+this.polarity_+"; RT-group: "+this.rtGroup_;
   }

  public String getTp()
  {
    return tp_;
  }

  public void setTp(String tp)
  {
    this.tp_ = tp;
  }

  public int getMaxFragments()
  {
    return maxFragments_;
  }

  public void setMaxFragments(int maxFragments)
  {
    this.maxFragments_ = maxFragments;
  }
  
  
  public int getMaxAllFrags()
  {
    return maxAllFrags_;
  }

  public void setMaxAllFrags(int maxAllFrags)
  {
    this.maxAllFrags_ = maxAllFrags;
  }

	public String getAmbiguousPartners()
	{
		return ambiguousPartners_;
	}
  
  
  
  
}
