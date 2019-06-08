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

import java.io.IOException;
import java.util.Hashtable;
import java.util.Set;
import java.util.Vector;

import at.tugraz.genome.lda.exception.NoRuleException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.msn.RulesContainer;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.QuantVO;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.quantification.CgProbe;

/**
 * object holding information about shared peaks - which are peaks that match to more than one analyte
 * the information about the peaks is stored in the SharedPeakContributionVOs
 * @author Juergen Hartler
 *
 */
public class SharedMS1PeakVO
{
  
  /** this vector contains the information about the shared identifications*/
  private Vector<SharedPeakContributionVO> contributions_;
  
  /** constructor initializing the required objects*/
  public SharedMS1PeakVO(){
    contributions_ = new Vector<SharedPeakContributionVO>();
  }
  
  /**
   * adds information about a contribution from a sharing analyte
   * @param quantVO the object containing quantitation instructions
   * @param set contains information about the MS1 peak and the MSn identification
   */
  public void addSharedInstance(QuantVO quantVO, LipidParameterSet set){
    contributions_.add(new SharedPeakContributionVO(quantVO,set));
  }
  
  /**
   * 
   * @return true if there are any MSn identifications present
   */
  public boolean areThereMS2Hits(){
    for (SharedPeakContributionVO contrib : contributions_){
      if (contrib.getSet() instanceof LipidomicsMSnSet)
        return true;
    }
    return false;
  }
  
  /**
   * 
   * @return the objects which contain information about the shared identifications 
   */
  public Vector<SharedPeakContributionVO> getPartners(){
    return contributions_;
  }
  
  /**
   * all of the added contributing MSn identifications are checked if they contain any 
   * distinct fragments (fragments that are not identified by any other partners)
   */
  public void checkForDistinctFragments(){
    for (int i=0; i!=contributions_.size();i++){
      SharedPeakContributionVO shared1 = contributions_.get(i);
      if (!(shared1.getSet() instanceof LipidomicsMSnSet)) continue;
      LipidomicsMSnSet set1 = (LipidomicsMSnSet) shared1.getSet();
      Hashtable<String,CgProbe> frags1 = buildFragmentHash(set1);
      Hashtable<String,CgProbe> distinctFragments = new Hashtable<String,CgProbe>();
      Hashtable<Integer,Boolean> msLevels = new Hashtable<Integer,Boolean>();
      for (String key1 : frags1.keySet()){
        CgProbe probe1 = frags1.get(key1);
        if (probe1.isFromOtherSpecies()) continue;
        msLevels.put(probe1.getMsLevel(),true);
        boolean distinct = true;
        for (int j=0; j!=contributions_.size(); j++){
          SharedPeakContributionVO shared2 = contributions_.get(j);
          if (i==j || !(shared2.getSet() instanceof LipidomicsMSnSet)) continue;
          LipidomicsMSnSet set2 = (LipidomicsMSnSet) shared2.getSet();
          Hashtable<String,CgProbe> frags2 = buildFragmentHash(set2);
          boolean foundOther = false;
          for (String key2 : frags2.keySet()){
            CgProbe probe2 = frags2.get(key2);
            if (probe2.isFromOtherSpecies()) continue;
            if (isInsideMz(probe1, probe2)){
              foundOther = true;
              break;
            }
          }
          if (foundOther){
            distinct = false;
            break;
          }
        }
        if (distinct) distinctFragments.put(key1, probe1);
      }
      // now check which fragments originate from other partners
      Hashtable<String,CgProbe> fragsFromOthers = new Hashtable<String,CgProbe>();
      for (int j=0; j!=contributions_.size(); j++){
        SharedPeakContributionVO shared2 = contributions_.get(j);
        if (i==j || !(shared2.getSet() instanceof LipidomicsMSnSet)) continue;
        LipidomicsMSnSet set2 = (LipidomicsMSnSet) shared2.getSet();
        Hashtable<String,CgProbe> frags2 = buildFragmentHash(set2);
        for (String key2 : frags2.keySet()){
          CgProbe probe2 = frags2.get(key2);
          boolean found = false;
          for (String key1 : frags1.keySet()){
            CgProbe probe1 = frags1.get(key1);
            if (isInsideMz(probe1, probe2)){
              found = true;
              break;
            }
          }
          if (found) continue;
          for (CgProbe probe1 : fragsFromOthers.values()){
            if (isInsideMz(probe1, probe2)){
              found = true;
              break;
            }
          }
          if (!found) fragsFromOthers.put(shared2.getQuantVO().getAnalyteClass()+shared2.getSet().getNameStringWithoutRt()+"x"+key2, probe2);
        }
      }
      if (distinctFragments.size()>0){
//        System.out.println(shared1.getQuantVO().getAnalyteClass()+"-"+shared1.getSet().getNameString()+";"+fragsFromOthers.size());
        shared1.setDistinctFragments(distinctFragments);
      }
      shared1.setFoundFragments(frags1);
      shared1.setFragsFromOthers(fragsFromOthers);
      shared1.setMsLevels(msLevels);
    }
  }
  
  /**
   * builds a hash table containing distinct fragments for one MSn identification (shared m/z values are counted only once)
   * @param set contains information about the MS1 peak and the MSn identification
   * @return the hash table containing distinct fragments for one MSn identification (shared m/z values are counted only once)
   */
  private Hashtable<String,CgProbe> buildFragmentHash(LipidomicsMSnSet set){
    Hashtable<String,CgProbe> frags = new Hashtable<String,CgProbe>();
    for (String fragName : set.getHeadGroupFragments().keySet()) frags.put(fragName, set.getHeadGroupFragments().get(fragName));
    Hashtable<String,Hashtable<String,CgProbe>> chainFragments = set.getChainFragments();
    for (String faName : chainFragments.keySet()){
      Hashtable<String,CgProbe> chainFrags = chainFragments.get(faName);
      for (String fragName : chainFrags.keySet()){
        CgProbe probe = chainFrags.get(fragName);
        String otherProbe = null;
        // check if there is already a probe with this m/z
        for (String other : frags.keySet()){
          if (isInsideMz(probe, frags.get(other))){
            otherProbe = other;
            break;
          }
        }
        if (otherProbe==null || probe.Area>frags.get(otherProbe).Area) {
          if (otherProbe!=null)
            frags.remove(otherProbe);
          frags.put(faName+";"+fragName, probe);
        }
      }
    }
    return frags;
  }
  
  /**
   * check if the m/z values of two different fragments overlap -> if so, they are not distinct
   * @param probe1 the first fragment
   * @param probe2 the second fragment
   * @return true if the fragments overlap -> if so, they are not distinct
   */
  private boolean isInsideMz(CgProbe probe1, CgProbe probe2){
    boolean inside = false;
    if (probe1.Mz<probe2.Mz){
      if (probe1.Mz+probe1.UpperMzBand>probe2.Mz-probe2.LowerMzBand)
        inside = true;
    } else {
      if (probe2.Mz+probe2.UpperMzBand>probe1.Mz-probe1.LowerMzBand)
        inside = true;
    }
    return inside;
  }
  
  /**
   * 
   * @return true if there is a contributing species that contains distinct fragments
   */
  public boolean hasAnyPartnerDistinctFragments(){
    for (SharedPeakContributionVO contr :  contributions_){
      if (contr.hasDistinctFragments()) return true;
    }
    return false;
  }
  
  /**
   * removes contributing species from the vector (usually the ones without distinct fragments)
   * @param partnersToRemove vector containing the positions of the partners to remove
   */
  public void removePartners(Vector<Integer> partnersToRemove){
    for (int i=(partnersToRemove.size()-1); i!=-1; i--){
      contributions_.remove(partnersToRemove.get(i).intValue());
    }
  }
  
  /**
   * checks if the spectral contribution of the distinct peaks of this species  is lower than a threshold -
   * compared to the peak contributing most 
   * @param threshold the relative cutoff threshold
   * @param ofInterest the species to be tested
   * @return true if the spectral contribution of the distinct peaks of this species  is lower than the threshold
   */
  public boolean isSpectrumContributionMuchLower (float threshold, SharedPeakContributionVO ofInterest){
    boolean isMuchLower = false;
    float totalContribution = 0f;
    for (CgProbe probe : ofInterest.getDistinctFragments().values()){
      totalContribution += probe.Area;
    }
    for (SharedPeakContributionVO other : contributions_){
      if (ofInterest.getQuantVO().equals(other.getQuantVO())) continue;
      float otherContrib = 0f;
      for (CgProbe probe : other.getDistinctFragments().values()){
        otherContrib += probe.Area;
      }
//      System.out.println(threshold+";"+totalContribution/(totalContribution+otherContrib));
      if (threshold>(totalContribution/(totalContribution+otherContrib))){
        isMuchLower = true;
        break;
      }
    }
//    System.out.println(isMuchLower);
    return isMuchLower;
  }
  
  /**
   * checks whether classes/adducts are present in the shared peak instances
   * @param adducts the adducts to be checked
   * @return true when classes/adducts are present
   */
  public boolean areAnyOfTheseAdductsPresent(Set<String> adducts) {
    boolean there = false;
    for (SharedPeakContributionVO contr :  contributions_){
      if (adducts.contains(StaticUtils.getRuleName(contr.getQuantVO().getAnalyteClass(), contr.getQuantVO().getModName()))) {
        there = true;
        break;
      }
    }
    return there;
  }
  
  /**
   * 
   * @return true if there is a contributing species that contains distinct fragments
   */
  public boolean haveAllChooseOnRtSetToTrue(){
    boolean allTrue = contributions_.size()>1;
    //for testing
    for (SharedPeakContributionVO contr :  contributions_){      
      try {
        if (!RulesContainer.choseMoreLikelyRtWhenEqualMSn(StaticUtils.getRuleName(contr.getQuantVO().getAnalyteClass(),contr.getQuantVO().getModName()))) {
          allTrue = false;
          break;
        }
      }
      catch (RulesException | NoRuleException | IOException
          | SpectrummillParserException e) {
        allTrue=false;
      }
    }
    return allTrue;
  }
}
