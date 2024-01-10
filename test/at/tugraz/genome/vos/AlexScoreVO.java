package at.tugraz.genome.vos;

public class AlexScoreVO
{
  
  private String lipidClass_;
  private String lipidSpecies_;
  private String molSpecies_;
  private short identificationType_;
  private double score_;
  private String adduct_;
  private boolean polarity_;
  private String rtGroup_;
  private double prophetScore_;
  private String tp_;
  private String comment_;
  
  public final static short IDENT_TYPE_SPECIES = 0;
  public final static short IDENT_TYPE_MOL_SPECIES = 1;
  
  
  public AlexScoreVO(String lipidClass, String lipidSpecies,
      String molSpecies, short identificationType, double score,
      String adduct, boolean polarity, String rtGroup, String comment)
  {
    super();
    this.lipidClass_ = lipidClass;
    this.lipidSpecies_ = lipidSpecies;
    this.molSpecies_ = molSpecies;
    this.identificationType_ = identificationType;
    this.score_ = score;
    this.adduct_ = adduct;
    this.polarity_ = polarity;
    this.rtGroup_ = rtGroup;
    this.prophetScore_ = 0.0d;
    this.tp_ = "";
    this.comment_ = comment;
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
     return "Score: "+this.score_+"; Class: "+this.lipidClass_+"; Species: "+this.lipidSpecies_+";"+this.molSpecies_+"; Adduct: "+this.adduct_+
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
  
  
  
}
