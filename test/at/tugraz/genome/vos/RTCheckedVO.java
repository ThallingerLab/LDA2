package at.tugraz.genome.vos;

public class RTCheckedVO
{

  private String adduct_;
  private String score_;
  private String fragType_;
  private String lClass_;
  private String species_;
  private String molSpec_;
  private String polarity_;
  private String rtGroup_;
  private String truePos_;
  private String comment_;
  
  
  
  public RTCheckedVO(String adduct, String score, String fragType,
      String lClass, String species, String molSpec, String polarity,
      String rtGroup, String truePos, String comment)
  {
    super();
    this.adduct_ = adduct;
    this.score_ = score;
    this.fragType_ = fragType;
    this.lClass_ = lClass;
    this.species_ = species;
    this.molSpec_ = molSpec;
    this.polarity_ = polarity;
    this.rtGroup_ = rtGroup;
    this.truePos_ = truePos;
    this.comment_ = comment;
  }



  public String getAdduct()
  {
    return adduct_;
  }



  public void setAdduct(String adduct)
  {
    this.adduct_ = adduct;
  }



  public String getScore()
  {
    return score_;
  }



  public void setScore(String score)
  {
    this.score_ = score;
  }



  public String getFragType()
  {
    return fragType_;
  }



  public void setFragType(String fragType)
  {
    this.fragType_ = fragType;
  }



  public String getlClass()
  {
    return lClass_;
  }



  public void setlClass_(String lClass)
  {
    this.lClass_ = lClass;
  }



  public String getSpecies()
  {
    return species_;
  }



  public void setSpecies(String species)
  {
    this.species_ = species;
  }



  public String getMolSpec()
  {
    return molSpec_;
  }



  public void setMolSpec(String molSpec)
  {
    this.molSpec_ = molSpec;
  }



  public String getPolarity()
  {
    return polarity_;
  }



  public void setPolarity_(String polarity)
  {
    this.polarity_ = polarity;
  }



  public String getRtGroup()
  {
    return rtGroup_;
  }



  public void setRtGroup_(String rtGroup)
  {
    this.rtGroup_ = rtGroup;
  }



  public String getTruePos()
  {
    return truePos_;
  }



  public void setTruePos_(String truePos)
  {
    this.truePos_ = truePos;
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
