/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp  
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

package at.tugraz.genome.lda.vos;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Set;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.utils.StaticUtils;
import javafx.util.Pair;

/**
 * Value object holding data required for the assignment of a double bond position
 * 
 * @author Leonida M. Lamp
 *
 */
public class DoubleBondPositionVO implements Comparable<DoubleBondPositionVO>
{
  /** Vector of FattyAcidVOs the molecular species is composed of */
  Vector<FattyAcidVO> chainCombination_;
  /** the expected retention time */
  float expectedRetentionTime_;
  /** the accuracy of this retention time match */
  int accuracy_;
  /** molecular species String without double bond positions */
  String molecularSpecies_;
  /** whether this double bond position is assigned */
  boolean isAssigned_;
  /** To optionally add Exp, RT pairs TODO: only for analysis, remove after */
  Hashtable<Pair<String,Integer>,Float> experimentRTLookup_ = new Hashtable<Pair<String,Integer>,Float>();
  
  /**
   * Constructor for a value object holding data required for the assignment of a double bond position
   * @param chainCombination Vector of FattyAcidVOs the molecular species is composed of
   * @param expectedRetentionTime the expected retention time
   */
  public DoubleBondPositionVO(
      Vector<FattyAcidVO> chainCombination, float expectedRetentionTime) 
  {
    this(chainCombination, expectedRetentionTime,0,null);
  }
  
  /**
   * Constructor for a value object holding data required for the assignment of a double bond position
   * @param chainCombination Vector of FattyAcidVOs the molecular species is composed of
   * @param expectedRetentionTime the expected retention time
   * @param accuracy the accuracy of this retention time match
   * @param molecularSpecies molecular species String without double bond positions
   */
  public DoubleBondPositionVO(
      Vector<FattyAcidVO> chainCombination, float expectedRetentionTime, int accuracy, String molecularSpecies) 
  {
    this(chainCombination, expectedRetentionTime, accuracy, molecularSpecies, false);
  }
  
  /**
   * Constructor for a value object holding data required for the assignment of a double bond position
   * @param chainCombination Vector of FattyAcidVOs the molecular species is composed of
   * @param expectedRetentionTime the expected retention time
   * @param accuracy the accuracy of this retention time match
   * @param molecularSpecies molecular species String without double bond positions
   * @param isAssigned whether this double bond position is assigned
   */
  public DoubleBondPositionVO(
      Vector<FattyAcidVO> chainCombination, float expectedRetentionTime, int accuracy, String molecularSpecies, boolean isAssigned) 
  {
    this.chainCombination_ = chainCombination;
    this.expectedRetentionTime_ = expectedRetentionTime;
    this.accuracy_ = accuracy;
    this.molecularSpecies_ = molecularSpecies;
    this.isAssigned_ = isAssigned;
  }
  
  /**
   * Copy Constructor
   * @param that DoubleBondPositionVO to be copied
   */
  public DoubleBondPositionVO(DoubleBondPositionVO that) 
  {
    this(getDeepCopyOfChainCombination(that.getChainCombination()), 
        that.getExpectedRetentionTime(), 
        that.getAccuracy(), 
        that.getMolecularSpecies(), 
        that.getIsAssigned());
  }
  
  /**
   * Creates a deep copy of a Vector of FattyAcidVOs
   * @param chainCombination the Vector of FattyAcidVOs
   * @return a deep copy of a Vector of FattyAcidVOs
   */
  private static Vector<FattyAcidVO> getDeepCopyOfChainCombination(Vector<FattyAcidVO> chainCombination)
  {
    Vector<FattyAcidVO> newChainCombination = new Vector<FattyAcidVO>();
    for (FattyAcidVO fattyAcid : chainCombination) {
      newChainCombination.add(new FattyAcidVO(fattyAcid));
    }
    return newChainCombination;
  }
  
  /**
   * Compares this object with the specified object for order. 
   * First, accuracy is compared (species with a higher accuracy are given priority)
   * Objects with the same accuracy are compared by their molecular species without potential prefixes or double bond assignments
   * Objects with the same molecular species are compared by their expected retention time.
   * Objects with the same expected retention time are then compared by the sum of their omega positions (Objects with a lower sum are given priority).
   */
  public int compareTo(DoubleBondPositionVO anotherVO) 
  {
  	return Comparator
  			.comparing(DoubleBondPositionVO::getAccuracy).reversed()
  			.thenComparing((DoubleBondPositionVO vo) -> vo.getEncodedDetailed(false, false))
  			.thenComparing(DoubleBondPositionVO::getExpectedRetentionTime)
  			.thenComparing(DoubleBondPositionVO::getOmegaSum)
  			.compare(this,anotherVO);
  }
  
  /**
   * @return the molecular species including double bond position information at the position level as a human readable String
   */
  public String getDoubleBondPositionsHumanReadable() 
  {
    return getDoubleBondPositionsHumanReadable(LipidomicsConstants.EXPORT_ANALYTE_TYPE_POSITION);
  }
  
  /**
   * @param snLevel 		only when this parameter equals at least LipidomicsConstants.EXPORT_ANALYTE_TYPE_POSITION, 
   * 										available sn position information is included in the molecular species String
   * @return the molecular species including double bond position information as a human readable String
   */
  public String getDoubleBondPositionsHumanReadable(int snLevel) 
  {
    boolean chainPositionsFixed = false;
    if (snLevel == LipidomicsConstants.EXPORT_ANALYTE_TYPE_POSITION) 
    {
      chainPositionsFixed = this.areChainPositionsFixed();
    }
    String doubleBondPositionsHumanReadable = null;
    try 
    {
      doubleBondPositionsHumanReadable = StaticUtils.getHumanReadableCombiName(this.getChainCombination(),
          Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding(),chainPositionsFixed);
    }
    catch (LipidCombinameEncodingException lcx) 
    {
      lcx.printStackTrace();
    }
    return doubleBondPositionsHumanReadable;
  }
  
  /**
   * Checks whether the molecular species String, if given, has known chain positions. 
   * If so, a check is performed whether this can also be applied when double bond position information is given.
   * @return whether chain positions are assigned for this molecular species with double bond position information
   */
  public boolean areChainPositionsFixed() 
  {
    if ((molecularSpecies_ != null) && molecularSpecies_.contains(LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS)) 
    {
      String[] individualChains = molecularSpecies_.split(LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS);
      if (individualChains.length > 1) 
      {
        if (areDoubleBondPositionsAssignedForAllChains()) 
        {
          return true;
        }
        Set<String> set = new HashSet<String>(Arrays.asList(individualChains));
        Object[] uniqueChains = set.toArray();
        if (individualChains.length == uniqueChains.length) 
        {
          return true;
        }
      }
    }
    return false;
  }
  
  /**
   * @return whether all chains this molecular species is composed of have assigned double bond positions
   */
  public boolean areDoubleBondPositionsAssignedForAllChains() 
  {
    for (FattyAcidVO fattyAcid : chainCombination_) 
    {
      if (fattyAcid.getOmegaPosition() < 0) 
      {
        return false;
      }
    }
    return true;
  }
  
  /**
   * @return the double bond position assignments for the chains in the order they are sorted in
   */
  public Vector<Integer> getPositionAssignmentPattern()
  {
  	Vector<Integer> pattern = new Vector<Integer>();
  	for (FattyAcidVO fattyAcid : chainCombination_) 
    {
  		pattern.add(fattyAcid.getOmegaPosition());
    }
  	return pattern;
  }
  
  private int getOmegaSum()
  {
  	Vector<Integer> pattern = this.getPositionAssignmentPattern();
  	int sum = 0;
  	for (int omega : pattern)
  	{
  		sum += omega;
  	}
  	return sum;
  }
  
  /**
   * @return the pattern in which the chains (in the order they are sorted in) are allowed to be permuted
   */
  public Vector<Boolean> getPermutationPattern()
  {
  	Vector<Boolean> pattern = new Vector<Boolean>();
  	for (int i=0;i<chainCombination_.size();i++)
  	{
  		FattyAcidVO fattyAcid1 = chainCombination_.get(i);
  		pattern.add(i,false);
  		for (FattyAcidVO fattyAcid2 : chainCombination_) 
      {
  			if (fattyAcid1 != fattyAcid2) //we do not want to compare identical objects
  			{
  				if (fattyAcid1.equalsNotConsideringOmegaPosition(fattyAcid2))
  				{
  					pattern.setElementAt(true, i);
  				}
  			}
      }
  	}
  	return pattern;
  }
  
  public String getEncoded()
  {
  	return StaticUtils.encodeLipidCombi(chainCombination_);
  }
  
  /**
   * @param includePrefix						true if potential prefixes should be included in the encoded string
   * @param includeOmegaPosition		true if potential omega positions should be included in the encoded string
   * @return
   */
  public String getEncodedDetailed(boolean includePrefix, boolean includeOmegaPosition)
  {
  	Vector<String> chainIds = new Vector<String>();
    for (FattyAcidVO chain : chainCombination_) chainIds.add(chain.getChainIdDetailed(includePrefix, includeOmegaPosition));
    return StaticUtils.encodeLipidCombiFromIds(chainIds);
  }
  
  public int getNumberOfCarbons()
  {
  	int num = 0;
  	for (FattyAcidVO fa : chainCombination_)
  	{
  		num += fa.getcAtoms();
  	}
  	return num;
  }
  
  public int getNumberOfDoubleBonds()
  {
  	int num = 0;
  	for (FattyAcidVO fa : chainCombination_)
  	{
  		num += fa.getDoubleBonds();
  	}
  	return num;
  }
  
  public int getNumberOfOH()
  {
  	int num = 0;
  	for (FattyAcidVO fa : chainCombination_)
  	{
  		num += fa.getOhNumber();
  	}
  	return num == 0 ? -1 : num;
  }
  
  
  /**
   * @param isAssigned true if this double bond position information should be assigned, false otherwise 
   */
  public void setIsAssigned(boolean isAssigned) 
  {
    this.isAssigned_ = isAssigned;
  }
  
  /**
   * @return whether this double bond position information is assigned
   */
  public boolean getIsAssigned() 
  {
    return this.isAssigned_;
  }
  
  /**
   * @param molecularSpecies the molecular species String without double bond positions
   */
  public void setMolecularSpecies(String molecularSpecies) 
  {
    this.molecularSpecies_ = molecularSpecies;
  }

  /**
   * @return the molecular species String without double bond positions
   */
  public String getMolecularSpecies() 
  {
    return this.molecularSpecies_;
  }
  
  /**
   * @param expectedRetentionTime the expected retention time in minutes
   */
  public void setExpectedRetentionTime(float expectedRetentionTime) 
  {
    this.expectedRetentionTime_ = expectedRetentionTime;
  }
  
  /**
   * @return the expected retention time in minutes
   */
  public float getExpectedRetentionTime() 
  {
    return this.expectedRetentionTime_;
  }
  
  /**
   * @param chainCombination the Vector of FattyAcidVOs the molecular species is composed of
   */
  public void setChainCombination(Vector<FattyAcidVO> chainCombination) 
  {
    this.chainCombination_ = chainCombination;
  }
  
  /**
   * @return the Vector of FattyAcidVOs the molecular species is composed of
   */
  public Vector<FattyAcidVO> getChainCombination()
  {
    return this.chainCombination_;
  }
  
  /**
   * @param accuracy the accuracy of this retention time match
   */
  public void setAccuracy(int accuracy) 
  {
    this.accuracy_ = accuracy;
  }
  
  /**
   * @return the accuracy of this retention time match
   */
  public int getAccuracy() 
  {
    return this.accuracy_;
  }

  /**
   * TODO just for analysis, remove after
   * @return
   */
	public Hashtable<Pair<String,Integer>,Float> getExperimentRTLookup()
	{
		return experimentRTLookup_;
	}

	/**
	 * TODO just for analysis, remove after
	 * @param experiment
	 * @param rt
	 */
	public void addToExperimentRTLookup(Pair<String,Integer> experiment, Float rt)
	{
		this.experimentRTLookup_.put(experiment, rt);
	}
  
}
