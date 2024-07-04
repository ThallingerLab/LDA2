/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2024 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER. 
 *  
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * by the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version.
 *  
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details. 
 *  
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * Please contact lda@genome.tugraz.at if you need additional information or 
 * have any questions.
 */

package at.tugraz.genome.vos;

import java.util.Comparator;
import java.util.Objects;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class RTCheckedVO implements Comparable<RTCheckedVO>
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
  
  
  public RTCheckedVO(RTCheckedVO other)
	{
  	this(other.getAdduct(), other.getScore(), other.getFragType(),
        other.getlClass(), other.getSpecies(), other.getMolSpec(), other.getPolarity(),
        other.getRtGroup(), other.getTruePos(), other.getComment());
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
  
  
  
  private Double getScoreDouble()
  {
  	return Double.parseDouble(score_);
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


	@Override
	public int hashCode()
	{
		return Objects.hash(adduct_, comment_, fragType_, lClass_, molSpec_,
				polarity_, rtGroup_, score_, species_, truePos_);
	}


	@Override
	public boolean equals(Object obj)
	{
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		RTCheckedVO other = (RTCheckedVO) obj;
		return Objects.equals(adduct_, other.adduct_)
				&& Objects.equals(comment_, other.comment_)
				&& Objects.equals(fragType_, other.fragType_)
				&& Objects.equals(lClass_, other.lClass_)
				&& Objects.equals(molSpec_, other.molSpec_)
				&& Objects.equals(polarity_, other.polarity_)
				&& Objects.equals(rtGroup_, other.rtGroup_)
				&& Objects.equals(score_, other.score_)
				&& Objects.equals(species_, other.species_)
				&& Objects.equals(truePos_, other.truePos_);
	}


	@Override
	public int compareTo(RTCheckedVO o)
	{
		return Comparator
  			.comparing(RTCheckedVO::getScoreDouble).reversed()
  			.thenComparing(RTCheckedVO::getlClass)
  			.thenComparing(RTCheckedVO::getSpecies)
  			.thenComparing(RTCheckedVO::getAdduct)
  			.thenComparing(RTCheckedVO::getMolSpec)
  			.thenComparing(RTCheckedVO::getFragType)
  			.thenComparing(RTCheckedVO::getRtGroup)
  			.thenComparing(RTCheckedVO::getTruePos)
  			.compare(this,o);
	}
  
  
}
