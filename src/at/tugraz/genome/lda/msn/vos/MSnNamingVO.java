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

import java.util.HashSet;
import java.util.Set;
import java.util.Vector;

import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.swing.Range;
import at.tugraz.genome.lda.utils.StaticUtils;

/**
 * value object holding the various MSn identification names and providing
 * methods for naming comparision based on the various levels of evidence
 * furthermore, the method retains the retention time range the identification was observed
 * this class is used for merging peaks
 * @author Juergen Hartler
 *
 */
public class MSnNamingVO
{
  
  /** the MSn evidence status*/
  private int status_;
  /** the MS1 identification name*/
  private String ms1Name_;
  /** the permuted options for FA names*/
  private Set<String> faNames_;
  /** the permuted options for position names (several if not all positions can be assigned)*/
  private Set<String> positionNames_;
  /** the peak RTs where this peak was observed*/
  private Set<String> rts_;
  /** a retention time limit where another higher level contradicting element was founnd - otherwise everything would be merged on the MS1 level*/
  private float higherLevelStop_;
  /** a retention time limit where another higher level contradicting element was founnd - otherwise everything would be merged on the MS1 level*/  
  private float higherLevelStart_;
  /** the range where it is possible to merge this identification to another one - depends on the merging time*/
  private Range mergingRange_;
  /** max time distanc to allow for a union of peaks*/
  private float mergingTime_;
  
  /**
   * constructor using an detected identification and the maximum amount of time for merging in minutes
   * @param set the MSn identification
   * @param mergingTime maximum amount of time for merging in minutes
   */
  public MSnNamingVO(LipidomicsMSnSet set, float mergingTime){
    this(set.getNameStringWithoutRt(),set.getMSnIdentificationNames(), set.getStatus(), set.getRt(), mergingTime);
  }
  
  /**
   * constructor using the ms1Name, the MSn names, the MSn evidence status, the RT of the identification, and the maximum amount of time for merging in minutes
   * @param ms1Name common name of the MS1 identification
   * @param nameObject object holding the MSn identification names
   * @param status MSn evidence status
   * @param rt retention time of the identification
   * @param mergingTime maximum amount of time for merging in minutes
   */
  @SuppressWarnings("unchecked")
  public MSnNamingVO(String ms1Name, Vector<Object> nameObject, int status, String rt, float mergingTime){
    this.status_ = status;
    this.mergingTime_ = mergingTime;
    higherLevelStop_ = Float.MAX_VALUE;
    higherLevelStart_ = 0f;
    rts_ = new HashSet<String>();
    this.addPermanentRt(rt);
    this.ms1Name_ = ms1Name;
    faNames_ = new HashSet<String>();
    positionNames_ = new HashSet<String>();
    if (status_>LipidomicsMSnSet.HEAD_GROUP_DETECTED){
      for (Object names : nameObject){
        if (names instanceof Vector){
          for (String name : (Vector<String>)names){
            positionNames_.add(name);
            detectAllPotentialFaCombinations(name);
          }
        }else{
          String name = (String)names;
          if (status_==LipidomicsMSnSet.POSITION_DETECTED && name.indexOf("/")!=-1) positionNames_.add(name);
          detectAllPotentialFaCombinations(name);
        }
      }
    }
  }
  
  /**
   * permutes fatty acids to all potential positions and constructs names out of this
   * @param name the currently used combined name for the various fatty acids
   */
  private void detectAllPotentialFaCombinations(String name){
    String nm = name.replaceAll("/", "_");
    String[] faArray = nm.split("_");
    Vector<String> fas = new Vector<String>();
    for (String fa : faArray) fas.add(fa);
    faNames_.addAll(StaticUtils.getPermutedChainNames(fas));
  }
  
  /**
   * adds another found retention time to the list and calculates a new merging range
   * @param rt the retention time of the additional hit
   */
  private void addPermanentRt(String rt){
    this.rts_.add(rt);
    calculateInfluenceRange();
  }
  
  /**
   * calculates a new potential merging range
   */
  private void calculateInfluenceRange(){
    float lowestRt = Float.MAX_VALUE;
    float highestRt = 0f;
    for (String rtString :rts_){
      float rTime = Float.parseFloat(rtString);
      if (rTime<lowestRt) lowestRt = rTime;
      if (rTime>highestRt) highestRt = rTime;
    }
    float start = lowestRt-mergingTime_;
    if (higherLevelStart_>start) start = higherLevelStart_;
    float stop = highestRt+mergingTime_;
    if (higherLevelStop_<stop) stop = higherLevelStop_;
    mergingRange_ = new Range(start,stop);
  }

  /**
   * 
   * @return MSn evidence status of this hit
   */
  public int getStatus()
  {
    return status_;
  }
  
  /**
   * checks if the retention time of another hit is inside the potential merging range
   * @param other the VO
   * @return true if the retention time of another hit is inside the potential merging range
   */
  public boolean insideRange(MSnNamingVO other){
    boolean insideRange = false;
    for (String rtString : other.rts_){
      if (insideRange(Float.parseFloat(rtString))){
        insideRange = true;
        break;
      }
    }
    return insideRange;
  }
  
  /**
   * checks if a retention time is inside the potential merging range
   * @param rt
   * @return true if a retention time is inside the potential merging range
   */
  public boolean insideRange(float rt){
    return this.mergingRange_.insideRange(rt);
  }
  
  /**
   * contains the identification at this certain MSn evidence level the same evidence, or does it contradict
   * @param status MSn evidence status
   * @param other the other MSn identification
   * @return true if the identification at this certain MSn evidence level the same evidence
   */
  public boolean hasSameEvidence(int status, MSnNamingVO other){
    if (status==LipidomicsMSnSet.POSITION_DETECTED){
      for (String name : positionNames_){
        for (String otherName : other.positionNames_){
          if (name.equalsIgnoreCase(otherName)){
            return true;
          }
        }
      }
    } else if (status==LipidomicsMSnSet.FRAGMENTS_DETECTED){
      for (String name : faNames_){
        for (String otherName : other.faNames_){
          if (name.equalsIgnoreCase(otherName))
            return true;
        }
      }
    } else if (status==LipidomicsMSnSet.HEAD_GROUP_DETECTED){
      if (this.ms1Name_.equalsIgnoreCase(other.ms1Name_)) return true;
    }
    return false;
  }
  
  /**
   * merges the evidence of another merging VO to this one
   * @param other the other VO
   */
  public void mergeWithOther(MSnNamingVO other){
   Range thisRange = getRtHardRange();
   Range otherRange = other.getRtHardRange();
   if (thisRange.getStart()<otherRange.getStart()){
     if (other.higherLevelStop_<this.higherLevelStop_) this.higherLevelStop_ = other.higherLevelStop_;
   } else {
     if (other.higherLevelStart_>this.higherLevelStart_) this.higherLevelStart_= other.higherLevelStart_;
   }
   for (String rt : other.rts_) this.addPermanentRt(rt);
   positionNames_.addAll(other.positionNames_);
   faNames_.addAll(other.faNames_);
 }
 
  /**
   * is a VO inside the really found RTs (not in the possible merging range)
   * @param other the other VO of this kind
   * @return true if another VO is inside the really found RTs
   */
 public boolean insideCombinedRtRange(MSnNamingVO other){
   if (this.rts_.size()<2) return false;
   Range range = getRtHardRange();
   for (String rtString : other.rts_){
     float rt = Float.parseFloat(rtString);
     if (range.insideRange(rt)) return true;
   }
   return false;  
 }
 
 /**
  * calculates the distance of this rt to the hard covered range
  * @param rt retention time
  * @return absolute distance to the hard range - returns 0 if inside range
  */
 public float absoluteDistanceToHardRange(float rt){
   Range range = getRtHardRange();
   if (range.insideRange(rt)) return 0f;
   if (rt<range.getStart()){
     return range.getStart()-rt;
   }else{
     return rt-range.getStop();
   }
 }
 
 /**
  * calculates the RT-Range that is covered by found retention times
  * @return the RT-Range that is covered by found retention times
  */
 private Range getRtHardRange(){
   float start = Float.MAX_VALUE;
   float stop = 0f;
   for (String rtString : rts_){
     float rt = Float.parseFloat(rtString);
     if (rt<start) start = rt;
     if (rt>stop) stop = rt;
   }
   Range range = new Range(start,stop);
   return range;
 }
 
 /**
  * if there is contradicting evidence, retention time limits are set that cannot be overriden by lower evidence merging processes
  * @param other contradictin detection
  */
 public void setHigherLevelLimits(MSnNamingVO other){
   Range thisRange = getRtHardRange();
   Range otherRange = other.getRtHardRange();
   if (thisRange.getStart()<otherRange.getStart()){
     float border = (thisRange.getStop()+otherRange.getStart())/2;
     this.higherLevelStop_ =  border;
     other.higherLevelStart_ = border;
     this.calculateInfluenceRange();
     other.calculateInfluenceRange();
   }else{
     float border = (otherRange.getStop()+thisRange.getStart())/2;
     this.higherLevelStart_ = border;
     other.higherLevelStop_ = border;
     this.calculateInfluenceRange();
     other.calculateInfluenceRange();
   }
 }

 /**
  * 
  * @return MS1 identification name
  */
  public String getMs1Name(){
    return ms1Name_;
  }
 
  public String toString(){
    String result = "Range: "+this.mergingRange_.getStart()+"-"+this.mergingRange_.getStop()+"\n";
    for (String rt : this.rts_){
      result+=rt+"\n";
    }
    return result;
  }
  
}
