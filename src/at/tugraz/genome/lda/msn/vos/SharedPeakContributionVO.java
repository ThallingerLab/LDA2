/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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

package at.tugraz.genome.lda.msn.vos;

import java.util.Hashtable;

import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.vos.QuantVO;
import at.tugraz.genome.maspectras.quantification.CgProbe;

/**
 * Contains information about a sharing species of an MS1 peak identification
 * @author Juergen Hartler
 *
 */
public class SharedPeakContributionVO
{
  /** a hash table containing the distinct fragments of an identification*/
  private Hashtable<String,CgProbe> distinctFragments_;
  
  /** a hash table containing the fragments that originate from other species*/
  private Hashtable<String,CgProbe> fragsFromOthers_;
  
  /** a hash table containing all found fragments - for calculation of spectrum coverage*/
  private Hashtable<String,CgProbe> foundFragments_;
  
  /** the MS levels the found fragments require - for calculation of spectrum coverage*/
  private Hashtable<Integer,Boolean> msLevels_;
  
  /** the object containing quantitation instructions*/
  private QuantVO quantVO_;
  
  /** contains information about the MS1 peak and the MSn identification*/
  private LipidParameterSet set_;
  
  /** if a peak split has to be removed because a split partner has a wrong retention time, the unsplit peak version is stored*/
  private LipidParameterSet oldSet_;

  /**
   * constructor storing the minimum information required for a sharing species
   * @param quantVO the object containing quantitation instructions
   * @param set contains information about the MS1 peak and the MSn identification
   */
  public SharedPeakContributionVO(QuantVO quantVO, LipidParameterSet set)
  {
    this.quantVO_ = quantVO;
    this.set_ = set;
    this.distinctFragments_  = new Hashtable<String,CgProbe>();
    oldSet_ = null;
  }

  /**
   * 
   * @return true if there are any distinct fragments detected (have to be set before by setDistinctFragments)
   */
  public boolean hasDistinctFragments()
  {
    return distinctFragments_.size()>0;
  }

  /**
   * 
   * @return the object containing quantitation instructions
   */
  public QuantVO getQuantVO()
  {
    return quantVO_;
  }

  /**
   * 
   * @return object that contains information about the MS1 peak and the MSn identification
   */
  public LipidParameterSet getSet()
  {
    return set_;
  }
  
  /**
   * 
   * @param distinctFragments sets the distinct fragments for this sharing species
   */
  public void setDistinctFragments(Hashtable<String,CgProbe> distinctFragments){
    this.distinctFragments_ = distinctFragments;
  }

  /**
   * 
   * @return a hash table containing the distinct fragments of an identification
   */
  public Hashtable<String,CgProbe> getDistinctFragments()
  {
    return distinctFragments_;
  }

  /**
   * 
   * @param set object that contains information about the MS1 peak and the MSn identification
   */
  public void setSet(LipidParameterSet set)
  {
    this.set_ = set;
  }

  /**
   * 
   * @return a hash table containing the fragments that originate from other species
   */
  public Hashtable<String,CgProbe> getFragsFromOthers()
  {
    return fragsFromOthers_;
  }

  /**
   * 
   * @param fragsFromOthers a hash table containing the fragments that originate from other species
   */
  public void setFragsFromOthers (Hashtable<String,CgProbe> fragsFromOthers)
  {
    this.fragsFromOthers_ = fragsFromOthers;
  }

  /**
   * 
   * @return the MS levels the found fragments require - for calculation of spectrum coverage
   */
  public Hashtable<Integer,Boolean> getMsLevels()
  {
    return msLevels_;
  }

  /**
   * 
   * @param msLevels the MS levels the found fragments require - for calculation of spectrum coverage
   */
  public void setMsLevels(Hashtable<Integer,Boolean> msLevels)
  {
    this.msLevels_ = msLevels;
  }

  /**
   * 
   * @return a hash table containing all found fragments - for calculation of spectrum coverage
   */
  public Hashtable<String,CgProbe> getFoundFragments()
  {
    return foundFragments_;
  }

  /**
   * 
   * @param foundFragments a hash table containing all found fragments - for calculation of spectrum coverage
   */
  public void setFoundFragments(Hashtable<String,CgProbe> foundFragments)
  {
    this.foundFragments_ = foundFragments;
  }
  
  /**
   * if a peak split has to be removed because a split partner has a wrong retention time the unsplit peak version is stored,
   * this method creates the copy of the unsplit version
   */
  public void makeCopyOfOldIdentification(){
    if (set_ instanceof LipidomicsMSnSet){
      this.oldSet_ = new LipidomicsMSnSet((LipidomicsMSnSet)set_);
    }else{
      this.oldSet_ = new LipidParameterSet(set_);
    }
  }

  /**
   * if a peak split has to be removed because a split partner has a wrong retention time the unsplit peak version is stored
   * @return the unsplit peak version
   */
  public LipidParameterSet getOldSet()
  {
    return oldSet_;
  }
  
  
}
