/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2019 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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

package at.tugraz.genome.lda.msn;

import java.io.IOException;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Set;
import java.util.Vector;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.HydroxylationEncodingException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.exception.NoRuleException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.msn.vos.SharedMS1PeakVO;
import at.tugraz.genome.lda.msn.vos.SharedPeakContributionVO;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.LipidomicsAnalyzer;
import at.tugraz.genome.lda.quantification.LipidomicsChromatogram;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.QuantVO;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.quantification.CgException;
import at.tugraz.genome.maspectras.quantification.CgProbe;

/**
 * class for disentangling peaks that originate from different lipid classes
 * @author Juergen Hartler
 *
 */
public class MSnPeakSeparator
{
  /** hash table containing the split peaks*/
  private Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> result_;
  /** object that holds MS data and can quantify fragments of interest*/
  private LipidomicsAnalyzer analyzer_;
  /** MS level of the quantitation*/
  private int msLevel_;
  /** if a peak split has to be removed because a split partner has a wrong retention time, the unsplit peak version is stored*/
  private Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> peaksBeforeSplit_;
  /** adducts where other adducts are required - have to be excluded from separating*/
  private Set<String> adductsThatRequireOtherAdduct_;

  
  /**
   * constructor containing the required fields
   * @param hitsAccordingToQuant hash table containing the peaks to split
   * @param analyzer object that holds MS data and can quantify fragments of interest
   * @param msLevel MS level of the quantitation errors from the quantitation process
   * @param adductsThatRequireOtherAdduct class/adduct combinations that are in this list exclude a split
   */
  public MSnPeakSeparator(Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> hitsAccordingToQuant, 
      Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> peaksBeforeSplit, LipidomicsAnalyzer analyzer, int msLevel,
      Set<String> adductsThatRequireOtherAdduct) {
    this.result_ = new Hashtable<QuantVO,Hashtable<String,LipidParameterSet>>();;
    for (QuantVO quant : hitsAccordingToQuant.keySet()) {
      this.result_.put(quant, new Hashtable<String,LipidParameterSet>(hitsAccordingToQuant.get(quant)));
    }
    this.analyzer_ = analyzer;
    this.msLevel_ = msLevel;
    peaksBeforeSplit_ = peaksBeforeSplit;
    adductsThatRequireOtherAdduct_ = adductsThatRequireOtherAdduct;
  }
  
  /**
   * this method splits the shared MS1 peaks according to the possibilities delivered by MS2 hits
   * @param hitsAccordingToQuant the results of the quantitation
   * @param analyzer the analyzer for quantifying hits
   * @param msLevel the level where the peak is oberlapping
   * @return the disentangled results
   * @throws CgException exception that is thrown when there is something wrong with the fragment detection
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  @SuppressWarnings("unchecked")
  public Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> disentagleSharedMS1Peaks() throws CgException, LipidCombinameEncodingException{
    //first detect which peaks are shared by different lipid classes
    Vector<SharedMS1PeakVO> sharedPeaks = detectSharedMS1PeakInstances(result_);
    Vector<Integer> sharedToRemove = new Vector<Integer>();
    //remove the ones where other adducts are required -> do not disentangle these ones
    for (int i=0; i!=sharedPeaks.size();i++){
      SharedMS1PeakVO shared = sharedPeaks.get(i);
      if (shared.areAnyOfTheseAdductsPresent(adductsThatRequireOtherAdduct_)) {
        sharedToRemove.add(i);
      }
    }
    for (int i=(sharedToRemove.size()-1);i!=-1;i--){
      sharedPeaks.remove(sharedToRemove.get(i).intValue());
    }
    
    sharedToRemove = new Vector<Integer>();
//    System.out.println(sharedPeaks.size());
    for (int i=0; i!=sharedPeaks.size();i++){
      SharedMS1PeakVO shared = sharedPeaks.get(i);
      shared.checkForDistinctFragments();
      if (shared.hasAnyPartnerDistinctFragments()){
        removeDetectionsWhereDistinctFragmentsAreMissing(shared,result_);
        if (shared.getPartners().size()==1){
          recheckMS2HitWithoutIsobar(shared,result_,analyzer_);
          sharedToRemove.add(i);
        } else {
          int previousPartnerSize = shared.getPartners().size()+1;
          while (previousPartnerSize>shared.getPartners().size()){
            previousPartnerSize = shared.getPartners().size();
            shared.checkForDistinctFragments();
            removeDetectionsWhereSpectrumCoverageWithoutSharedIsToLow(shared,result_,analyzer_);
          }
          //System.out.println("!!!!!!!!!!!!!!!!!! "+shared.getPartners().size());
          
          if (shared.getPartners().size()==1){
            recheckMS2HitWithoutIsobar(shared,result_,analyzer_);
            sharedToRemove.add(i);
          }
        }
      }
      if (!shared.hasAnyPartnerDistinctFragments()) {
        boolean chooseOneByRt = shared.haveAllChooseOnRtSetToTrue();
        for (SharedPeakContributionVO contr : shared.getPartners()){
          contr.getSet().setPercentalSplit(100f);
          if (chooseOneByRt) {
            contr.getSet().setChoseMoreLikelyRtWhenEqualMSn(true);
          }
        }
        sharedToRemove.add(i);
        if (shared.getPartners().size()>0) checkIfAllPartnerFragmentsFulfillSpectrumCoverage(shared,result_,analyzer_);
      }
    }
    for (int i=(sharedToRemove.size()-1);i!=-1;i--){
      sharedPeaks.remove(sharedToRemove.get(i).intValue());
    }
//    System.out.println("1. "+sharedPeaks.size());
    // now, there is a list of shared peaks left, where each contributor has at least one distinct
    // fragment; now, analytes with very minor contriubtions are removed
    float relIsobarSpecCovExclusion = LipidomicsConstants.getMs2IsobarSCExclusionRatio();
    sharedToRemove = new Vector<Integer>();
    for (int i=0; i!=sharedPeaks.size();i++){
      SharedMS1PeakVO shared = sharedPeaks.get(i);
      Vector<SharedPeakContributionVO> contributions = new Vector<SharedPeakContributionVO>(shared.getPartners());
      for (int j=0; j!=contributions.size(); j++){
        SharedPeakContributionVO contr = contributions.get(j);
        float exclValue = relIsobarSpecCovExclusion;
        try{
          String cutoff = RulesContainer.getIsobarExclusionRatio(StaticUtils.getRuleName(contr.getQuantVO().getAnalyteClass(),contr.getQuantVO().getModName()));
          if (cutoff!=null) exclValue = Float.parseFloat(cutoff);
        } catch(NoRuleException nrx){
        } catch (RulesException | IOException | SpectrummillParserException e) {
          e.printStackTrace();
        }
        if (shared.isSpectrumContributionMuchLower(exclValue,contr)){
          removePartner(contr,shared,result_);
          if (shared.getPartners().size()==1){
            recheckMS2HitWithoutIsobar(shared,result_,analyzer_);
            sharedToRemove.add(i);
          }
        }
      }
    }
    Hashtable<Integer,Integer> removed = new  Hashtable<Integer,Integer>();
    for (int i=(sharedToRemove.size()-1);i!=-1;i--){
      int index = sharedToRemove.get(i);
      if (removed.containsKey(index)) continue;
      sharedPeaks.remove(index);
      removed.put(index, index);
    }
    
//    System.out.println("2. "+sharedPeaks.size());

    // now, there is a list of shared peaks left, where each contributor has at least one distinct
    // fragment; now we try to figure out if there are peaks of the same lipid class, that are
    // several minutes away, and the spectrum coverage contribution of the distinct peak is quite
    // low compared to the partnering peak - then an accidental hit based on same fragments is assumed
    // and the peak is left to the other partner
    
    //for shotgun: only separation by distinct fragments is possible
    if (LipidomicsConstants.isShotgun()==LipidomicsConstants.SHOTGUN_TRUE){
      for (SharedMS1PeakVO shared : sharedPeaks){
        calculatePercentualSplitValueAccordingToMSn(shared);
        for (SharedPeakContributionVO contr : shared.getPartners()){
          result_.get(contr.getQuantVO()).put("", contr.getSet());
        }
      }
      return result_;
    }
    
    //first, the peaks are sorted by ascending retention time
    Vector<SharedMS1PeakVO> sortedSharedPeaks = new Vector<SharedMS1PeakVO>();
    Hashtable<Integer,Boolean> disentangledPeaks = new Hashtable<Integer,Boolean>();
    sharedToRemove = new Vector<Integer>();
    for (SharedMS1PeakVO shared : sharedPeaks){
      double rt1 = Double.parseDouble(shared.getPartners().get(0).getSet().getRt());
      int count=0;
      for (int i=0; i!=sortedSharedPeaks.size();i++){
        double rt2 = Double.parseDouble(sortedSharedPeaks.get(i).getPartners().get(0).getSet().getRt());
        if (rt1<rt2) break;
        count++;
      }
      if (count==sortedSharedPeaks.size()) sortedSharedPeaks.add(shared);
      else sortedSharedPeaks.add(count,shared);
      disentangledPeaks.put(sortedSharedPeaks.size()-1, false);
    }
    float minRtDifference = LipidomicsConstants.getMs2IsobaricOtherRtDifference();
    float relIsobarFarAreaSpecCovExclusion = LipidomicsConstants.getMs2IsobarSCFarExclusionRatio();
    // now I check the shared contributions, if there is a unique detection in a certain retention time direction
    Hashtable<QuantVO,Hashtable<Integer,Boolean>> checkedSets = new Hashtable<QuantVO,Hashtable<Integer,Boolean>>();
    for (int i=0; i!=sharedPeaks.size();i++){
      SharedMS1PeakVO shared = sharedPeaks.get(i);
//      String classes = "";
//      String rts = "";
//      for (SharedPeakContributionVO partner : shared.getPartners()){
//        classes += partner.getQuantVO().getAnalyteClass()+"_"+partner.getSet().getNameString()+"-";
//        rts += partner.getSet().getRt()+"/";
//      }
//      classes = classes.substring(0,classes.length()-1);
//      rts = rts.substring(0,rts.length()-1);
//      System.out.println("1. "+classes+" "+rts);
      Vector<SharedPeakContributionVO> contributions = new Vector<SharedPeakContributionVO>(shared.getPartners());
      for (int j=0; j!=contributions.size(); j++){
        SharedPeakContributionVO contr = contributions.get(j);
        float exclValue = relIsobarFarAreaSpecCovExclusion;
        float rtDiff = minRtDifference;
        try{
          String cutoff = RulesContainer.getIsobarFarExclusionRatio(StaticUtils.getRuleName(contr.getQuantVO().getAnalyteClass(),contr.getQuantVO().getModName()));
          if (cutoff!=null) exclValue = Float.parseFloat(cutoff);
          cutoff = RulesContainer.getIsobarFarRtDifference(StaticUtils.getRuleName(contr.getQuantVO().getAnalyteClass(),contr.getQuantVO().getModName()));
          if (cutoff!=null) rtDiff = Float.parseFloat(cutoff);
        } catch(NoRuleException nrx){
        } catch (RulesException | IOException | SpectrummillParserException e) {
          e.printStackTrace();
        }
        if (disentangledPeaks.containsKey(i) && disentangledPeaks.get(i)) continue;
        Hashtable<Integer,Boolean> checked = new Hashtable<Integer,Boolean>();
        if (checkedSets.containsKey(contr.getQuantVO())) checked = checkedSets.get(contr.getQuantVO());
        if (checked.containsKey(j)) continue;
        checked.put(j, true);
        checkedSets.put(contr.getQuantVO(), checked);
        if (!result_.containsKey(contr.getQuantVO())) continue;
        @SuppressWarnings("rawtypes")
        Vector uniqueAndShared = getUniqueAndSharedContributions(i,sharedPeaks,disentangledPeaks,contr,result_.get(contr.getQuantVO()));
        if (uniqueAndShared.size()<2 || ((Vector<LipidParameterSet>)uniqueAndShared.get(0)).size()==0) continue;
        Vector<LipidParameterSet> uniques = (Vector<LipidParameterSet>)uniqueAndShared.get(0);
        Vector<Integer> sharedInsts = (Vector<Integer>)uniqueAndShared.get(1);
        int direction = detectDirectionOfUniqueHits(Double.parseDouble(contr.getSet().getRt()), uniques);
        if (direction==0) continue;
        // if the unique hits have higher retention times, the other shared hits have to be checked first,
        // because otherwise no direction is assignable
        if (direction>0){
          Vector<Integer> otherSharedDesc = getSharedDescendingOrder(j,sharedInsts);
          for (Integer index : otherSharedDesc){
            if (disentangledPeaks.containsKey(index) && disentangledPeaks.get(index)) continue;
            SharedMS1PeakVO other = sharedPeaks.get(index);
            for (SharedPeakContributionVO otherContr : new Vector<SharedPeakContributionVO>(other.getPartners())){
              if (disentangledPeaks.containsKey(index) && disentangledPeaks.get(index)) continue;
              if (!otherContr.getQuantVO().equals(contr)) continue;
              if (checked.containsKey(index)) continue;
              checked.put(index, true);
              if (!result_.containsKey(otherContr.getQuantVO())) continue;
              @SuppressWarnings("rawtypes")
              Vector otherUniqueAndShared = getUniqueAndSharedContributions(index,sharedPeaks,disentangledPeaks,otherContr,result_.get(otherContr.getQuantVO()));
              if (otherUniqueAndShared.size()<2 || ((Vector<LipidParameterSet>)otherUniqueAndShared.get(0)).size()==0) continue;
              Vector<LipidParameterSet> otherUniques = (Vector<LipidParameterSet>)uniqueAndShared.get(0);
              int otherDirection = detectDirectionOfUniqueHits(Double.parseDouble(otherContr.getSet().getRt()), otherUniques);
              if (otherDirection!=direction) continue;
              if (isSufficentlyFarAway(rtDiff,otherContr,otherUniques) && other.isSpectrumContributionMuchLower(exclValue,otherContr)){
                removePartner(otherContr,other,result_);
                if (other.getPartners().size()==1){
                  recheckMS2HitWithoutIsobar(other,result_,analyzer_);
                  sharedToRemove.add(i);
                  disentangledPeaks.put(index, true);
                }
              }
            }
            
          }
          if (!result_.containsKey(contr.getQuantVO())) continue;
          uniques = (Vector<LipidParameterSet>)getUniqueAndSharedContributions(i,sharedPeaks,disentangledPeaks,contr,result_.get(contr.getQuantVO())).get(0);
        }
        if (isSufficentlyFarAway(rtDiff,contr,uniques) && shared.isSpectrumContributionMuchLower(exclValue,contr)){
          removePartner(contr,shared,result_);
          if (shared.getPartners().size()==1){
            recheckMS2HitWithoutIsobar(shared,result_,analyzer_);
            sharedToRemove.add(i);
            disentangledPeaks.put(i, true);
          }
        }
      }
    }
    Collections.sort(sharedToRemove);
    removed = new  Hashtable<Integer,Integer>();
    for (int i=(sharedToRemove.size()-1);i!=-1;i--){
      int index = sharedToRemove.get(i);
      if (removed.containsKey(index)) continue;
      sharedPeaks.remove(index);
      removed.put(index, index);
    }
//    System.out.println("3. "+sharedPeaks.size());
    
    //peaks that are still marked as shared, have to run through the separation routine
    float percentForBordersAbsolute = 5f;
    float percentForBordersRelative = 10f;
    
    LipidomicsChromatogram chrom = null;
    if (sharedPeaks.size()>0){
      QuantVO quant = sharedPeaks.iterator().next().getPartners().get(0).getQuantVO();
      chrom = analyzer_.readOneChromatogram((float)quant.getAnalyteMass(),msLevel_);
    }
    for (SharedMS1PeakVO shared : sharedPeaks){
      float start = Float.MAX_VALUE;
      float stop = 0f;
      Vector<LipidParameterSet> sets = new Vector<LipidParameterSet>();
      for (SharedPeakContributionVO contr : shared.getPartners()) sets.add(contr.getSet());
      float[] borders = analyzer_.calculatePercentualBorders(percentForBordersAbsolute,chrom,sets);
      start = borders[0];
      stop = borders[1];
      float[] bordersRelative = analyzer_.calculatePercentualBorders(percentForBordersRelative,chrom,sets);
      float startRelative = bordersRelative[0];
      float stopRelative = bordersRelative[1];
      
      String classes = "";
      String rts = "";
      for (SharedPeakContributionVO partner : shared.getPartners()){
        classes += partner.getQuantVO().getAnalyteClass()+"_"+partner.getSet().getNameString()+"-";
        rts += partner.getSet().getRt()+"/";
      }
      classes = classes.substring(0,classes.length()-1);
      rts = rts.substring(0,rts.length()-1);
      // here the peak sharing algorithm is applied
      if (shared.getPartners().size()==2 && analyzer_.countMSnSpectraOfRegion(start, stop,2)>1){
        Vector<SharedPeakContributionVO> contrs = MSnAnalyzer.splitTwoIsobaricPeaks(analyzer_,start,stop,startRelative,stopRelative,shared,msLevel_);
        float zeroIsoAreaOne = 0f;
        float zeroIsoAreaTwo = 0f;
        if (contrs!=null && contrs.size()==2){
          
          if (contrs.get(0).getSet()!=null && contrs.get(0).getSet().getIsotopicProbes()!=null && contrs.get(0).getSet().getIsotopicProbes().size()>0){
            for (CgProbe probe : contrs.get(0).getSet().getIsotopicProbes().get(0)){
              zeroIsoAreaOne += probe.Area;
            }
          }
          if (contrs.get(1).getSet()!=null && contrs.get(1).getSet().getIsotopicProbes()!=null && contrs.get(1).getSet().getIsotopicProbes().size()>0){
            for (CgProbe probe : contrs.get(1).getSet().getIsotopicProbes().get(0)){
              zeroIsoAreaTwo += probe.Area;
            }
          }
        }
        // peak was separatable
        if (contrs!=null && contrs.size()==2 && zeroIsoAreaOne>0f && zeroIsoAreaTwo>0f){
          if (contrs.get(0).getSet()!=null){
            result_.get(contrs.get(0).getQuantVO()).put(contrs.get(0).getSet().getRt(), contrs.get(0).getSet());
            peaksBeforeSplit_.get(contrs.get(0).getQuantVO()).put(contrs.get(0).getSet().getRt(), contrs.get(0).getOldSet());
          }
          if (contrs.get(1).getSet()!=null){
            result_.get(contrs.get(1).getQuantVO()).put(contrs.get(1).getSet().getRt(), contrs.get(1).getSet());
            peaksBeforeSplit_.get(contrs.get(1).getQuantVO()).put(contrs.get(1).getSet().getRt(), contrs.get(1).getOldSet());
          }
        // no MS1 separation possible -> split according to MS2 intensities
        } else {
          calculatePercentualSplitValueAccordingToMSn(shared);
          for (SharedPeakContributionVO contr : shared.getPartners()){
            result_.get(contr.getQuantVO()).put(contr.getSet().getRt(), contr.getSet());
          }

        }
      // here the peak intensity is split according to the distinct fragments TODO: has to be implemented
      } else {
        calculatePercentualSplitValueAccordingToMSn(shared);
//        System.out.println(classes+" "+rts);
        for (SharedPeakContributionVO contr : shared.getPartners()){
          result_.get(contr.getQuantVO()).put(contr.getSet().getRt(), contr.getSet());
        }
      }
    }
    return result_;
  }


  /**
   * 
   * @param hitsAccordingToQuant quantitation hits that are assigned to QuantVOs (first key); second key is the retention time
   * @return Vector of SharedMS1PeakVO which are peaks that match to more than one analyte, and contain SharedPeakContributionVO
   *         which hold the QuantVO, the LipidParameterSet, and information about distinct fragments
   */
  public static Vector<SharedMS1PeakVO> detectSharedMS1PeakInstances(Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> hitsAccordingToQuant){
    Vector<SharedMS1PeakVO> sharedPeaks = new Vector<SharedMS1PeakVO>();
    Vector<QuantVO> quants = new Vector<QuantVO>(hitsAccordingToQuant.keySet());
    Hashtable<Integer,Hashtable<String,String>> alreadyAdded = new Hashtable<Integer,Hashtable<String,String>>();
    for (int i1=0; i1!=quants.size(); i1++){
      QuantVO quant1 = quants.get(i1);
      Hashtable<String,LipidParameterSet> hits1 =  hitsAccordingToQuant.get(quant1);
      Hashtable<String,String> added = new Hashtable<String,String>();
      if (alreadyAdded.containsKey(i1)) added = alreadyAdded.get(i1);
      for (String rt1 : hits1.keySet()){
        if (added.containsKey(rt1)) continue;
        LipidParameterSet set1 = hits1.get(rt1);        
        SharedMS1PeakVO sharedPeak = new SharedMS1PeakVO();
        sharedPeak.addSharedInstance(quant1, set1);
        boolean otherAdded = false;
        for (int i2=i1+1; i2<quants.size(); i2++){
          Hashtable<String,String> added2 = new Hashtable<String,String>();
          if (alreadyAdded.containsKey(i2)) added2 = alreadyAdded.get(i2);
          QuantVO quant2 = quants.get(i2);
          Hashtable<String,LipidParameterSet> hits2 =  hitsAccordingToQuant.get(quant2);
          for (String rt2 : hits2.keySet()){
            if (added2.containsKey(rt2)) continue;
            LipidParameterSet set2 = hits2.get(rt2);
            if (rt1.equalsIgnoreCase(rt2) || LipidomicsAnalyzer.isPeakCenterTheSame(set1, set2)){
              sharedPeak.addSharedInstance(quant2, set2);
              otherAdded = true;
              added2.put(rt2, rt2);
            }
          }
          alreadyAdded.put(i2, added2);
        }
        if (otherAdded && sharedPeak.areThereMS2Hits()){
          sharedPeaks.add(sharedPeak);
          added.put(rt1, rt1);
          alreadyAdded.put(i1, added);
        }
      }
    }
    return sharedPeaks;
  }

  
  /**
   * removes all contributing identifications from the SharedMS1PeakVO and the final results (hitsAccordingToQuant) that have no distinct fragments
   * @param shared the SharedMS1PeakVO which is a peak shared by more than one analyte
   * @param hitsAccordingToQuant quantitation hits that are assigned to QuantVOs (first key); second key is the retention time - is the result
   */
  private void removeDetectionsWhereDistinctFragmentsAreMissing(SharedMS1PeakVO shared, Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> hitsAccordingToQuant){
    Vector<Integer> partnersToRemove = new Vector<Integer>();
    Vector<SharedPeakContributionVO> partners = shared.getPartners();
    for (int i=0; i!=partners.size(); i++){
      SharedPeakContributionVO partner = partners.get(i);
      if (!partner.hasDistinctFragments()){
        Hashtable<String,LipidParameterSet> sets = hitsAccordingToQuant.get(partner.getQuantVO());
        String rt = "";
        if (LipidomicsConstants.isShotgun()!=LipidomicsConstants.SHOTGUN_TRUE)
          rt = partner.getSet().getRt();
        sets.remove(rt);
        if (sets.size()==0) hitsAccordingToQuant.remove(partner.getQuantVO());
        partnersToRemove.add(i);
      }
    }
    shared.removePartners(partnersToRemove);
  }

 
  /**
   * removes a partner from set of contributing species
   * @param contr the partner to be removed
   * @param shared the shared peak object containing the partners
   * @param hitsAccordingToQuant quantitation hits that are assigned to QuantVOs (first key); second key is the retention time - is the result
   */
  private void removePartner(SharedPeakContributionVO contr, SharedMS1PeakVO shared, Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> hitsAccordingToQuant){
    Vector<Integer> partnersToRemove = new Vector<Integer>();
    Vector<SharedPeakContributionVO> partners = shared.getPartners();
    for (int i=0; i!=partners.size(); i++){
      SharedPeakContributionVO partner = partners.get(i);
      if (partner.getQuantVO().equals(contr.getQuantVO())){
        Hashtable<String,LipidParameterSet> sets = hitsAccordingToQuant.get(partner.getQuantVO());
        String rt = "";
        if (LipidomicsConstants.isShotgun()!=LipidomicsConstants.SHOTGUN_TRUE)
          rt = partner.getSet().getRt();
        sets.remove(rt);
        if (sets.size()==0) hitsAccordingToQuant.remove(partner.getQuantVO());
        partnersToRemove.add(i);
      }
    }
    shared.removePartners(partnersToRemove);
  }


  /**
   * makes again an MS2 analysis with full MSn rules (absolute intensities are checked) - this is the case, if one peak is identified to belong to one species only
   * @param shared the SharedMS1PeakVO which is a peak shared by more than one analyte (in this case, only one should be left)
   * @param hitsAccordingToQuant quantitation hits that are assigned to QuantVOs (first key); second key is the retention time - is the result
   * @param analyzer the analyzer for quantifying hits
   * @throws CgException exception that is thrown when there is something wrong with the fragment detection
   */
  private void recheckMS2HitWithoutIsobar(SharedMS1PeakVO shared, Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> hitsAccordingToQuant,
      LipidomicsAnalyzer analyzer) throws CgException{
    for (SharedPeakContributionVO contr : shared.getPartners()){
      QuantVO oneSet = contr.getQuantVO();
      Hashtable<String,LipidParameterSet> quantsOfMod = hitsAccordingToQuant.get(oneSet);
      try {
        String rt = "";
        if (LipidomicsConstants.isShotgun()!=LipidomicsConstants.SHOTGUN_TRUE)
          rt = contr.getSet().getRt();
        //TODO: this is set to true in the meantime - maybe play around with caching in the future to improve calculation time
        MSnAnalyzer msnAnalyzer = new MSnAnalyzer(oneSet.getAnalyteClass(),oneSet.getModName(),contr.getSet(),analyzer,oneSet,true,false);  
        if (msnAnalyzer.checkStatus()==LipidomicsMSnSet.DISCARD_HIT) quantsOfMod.remove(rt);
        else quantsOfMod.put(rt,msnAnalyzer.getResult());
      }
      catch (RulesException | HydroxylationEncodingException | ChemicalFormulaException | LipidCombinameEncodingException e) {
        e.printStackTrace();
      }
      catch (IOException e) {
        e.printStackTrace();
      }
      catch (SpectrummillParserException e) {
        e.printStackTrace();
      }

    }
  }

  
  /**
   * here, for all participating species of a peak, a spectrum coverage is calculated, that does not include
   * fragments identified by other isobaric species; contributions are discarded that are below the class specific coverage value 
   * @param shared the SharedMS1PeakVO which is a peak shared by more than one analyte (in this case, only one should be left)
   * @param hitsAccordingToQuant quantitation hits that are assigned to QuantVOs (first key); second key is the retention time - is the result
   * @param analyzer the analyzer for retrieving intensities from the MS data
   */
  private void removeDetectionsWhereSpectrumCoverageWithoutSharedIsToLow(SharedMS1PeakVO shared, Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> hitsAccordingToQuant,
      LipidomicsAnalyzer analyzer){
    Vector<Integer> partnersToRemove = new Vector<Integer>();
    for (int i=0; i!=shared.getPartners().size();i++){
      SharedPeakContributionVO contr = shared.getPartners().get(i);
      boolean remove = false;
      try {
        //TODO: this is set to true in the meantime - maybe play around with caching in the future to improve calculation time
        MSnAnalyzer.prepareCachedSpectra(analyzer,contr.getSet(),true);
        float bpCutoff = (float)RulesContainer.getBasePeakCutoff(StaticUtils.getRuleName(contr.getQuantVO().getAnalyteClass(), contr.getQuantVO().getModName()));
        float coverage = (float)RulesContainer.getSpectrumCoverageMin(StaticUtils.getRuleName(contr.getQuantVO().getAnalyteClass(), contr.getQuantVO().getModName()));
        if (!MSnAnalyzer.isSpectrumCovered(analyzer, contr.getSet(), contr.getMsLevels(), new Vector<CgProbe>(contr.getFoundFragments().values()),
            new Vector<CgProbe>(contr.getFragsFromOthers().values()),bpCutoff, coverage))
          remove = true;
      }
      catch (RulesException | NoRuleException | IOException
          | SpectrummillParserException | CgException e) {
        e.printStackTrace();
        remove = true;
      }
      if (remove){
        String rt = "";
        if (LipidomicsConstants.isShotgun()!=LipidomicsConstants.SHOTGUN_TRUE)
          rt = contr.getSet().getRt();
        hitsAccordingToQuant.get(contr.getQuantVO()).remove(rt);
        partnersToRemove.add(i);
      }
    }
    shared.removePartners(partnersToRemove);
  }
  
  
  /**
   * checks if all the sum of found fragments (from all contributing partners) fulfills the lowest possible spectrum coverage - otherwise, the hits are removed
   * @param shared the SharedMS1PeakVO which is a peak shared by more than one analyte (in this case, only one should be left)
   * @param hitsAccordingToQuant quantitation hits that are assigned to QuantVOs (first key); second key is the retention time - is the result
   * @param analyzer the analyzer for retrieving intensities from the MS data
   */
  private void checkIfAllPartnerFragmentsFulfillSpectrumCoverage(SharedMS1PeakVO shared, Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> hitsAccordingToQuant,
      LipidomicsAnalyzer analyzer) {
    try {
      Vector<CgProbe> allFragments = new Vector<CgProbe>();
      
      for (int i=0; i!=shared.getPartners().size();i++){
        SharedPeakContributionVO contr = shared.getPartners().get(i);
        if (contr.getSet() instanceof LipidomicsMSnSet){
          allFragments.addAll(contr.getFoundFragments().values());
          allFragments.addAll(contr.getFragsFromOthers().values());
          break;
        }
      }
      float minCutoff = 1f;
      float lowestCov = 1f;
      Hashtable<Integer,Boolean> msLevels = new Hashtable<Integer,Boolean>();
      for (int i=0; i!=shared.getPartners().size();i++){
        SharedPeakContributionVO contr = shared.getPartners().get(i);
        if (!(contr.getSet() instanceof LipidomicsMSnSet)) continue;
        String className = contr.getQuantVO().getAnalyteClass();
        String modName = contr.getQuantVO().getModName();
        float bpCutoff;
          bpCutoff = (float)RulesContainer.getBasePeakCutoff(StaticUtils.getRuleName(className, modName));
        float coverage = (float)RulesContainer.getSpectrumCoverageMin(StaticUtils.getRuleName(className, modName));
        if (bpCutoff<minCutoff) minCutoff = bpCutoff;
        if (coverage<lowestCov) lowestCov = coverage;
        for (Integer level : contr.getMsLevels().keySet()){
          if (contr.getMsLevels().get(level)) msLevels.put(level, true);
        }
      }
      SharedPeakContributionVO contr = null;
      for (SharedPeakContributionVO contr1 : shared.getPartners()){
        if (contr1.getSet() instanceof LipidomicsMSnSet){
          contr = contr1;
          break;
        }
      }
    //TODO: this is set to true in the meantime - maybe play around with caching in the future to improve calculation time
      MSnAnalyzer.prepareCachedSpectra(analyzer,contr.getSet(),true);
      if (!MSnAnalyzer.isSpectrumCovered(analyzer, contr.getSet(), msLevels, allFragments, new Vector<CgProbe>(),minCutoff,lowestCov)){
        for (int i=0; i!=shared.getPartners().size();i++){
          SharedPeakContributionVO cont = shared.getPartners().get(i);
          if (!(cont.getSet() instanceof LipidomicsMSnSet)) continue;
          String rt = "";
          if (LipidomicsConstants.isShotgun()!=LipidomicsConstants.SHOTGUN_TRUE)
            rt = cont.getSet().getRt();
          hitsAccordingToQuant.get(cont.getQuantVO()).remove(rt);
        }
      }

    }
    catch (RulesException | NoRuleException | IOException
        | SpectrummillParserException | CgException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }


  /**
   * returns two vectors - the first is of type LipidParameterSet and contains peaks that were uniquely identified (no shared instances) -
   * the second one is a list of indices that name other shared instances that overlap as well
   * @param currentSharedPeak the index of the currently shared peak under observation
   * @param sharedPeaks the vector containing all shared peak objects
   * @param disentangledPeaks boolean hash defining if shared peak instances were already entangled and only one species is left
   * @param contr the current contributing partner where the unique and shared instances have to be found
   * @param rtHits all found currently valid hits
   * @return
   */
  @SuppressWarnings({ "rawtypes", "unchecked" })
  private Vector getUniqueAndSharedContributions(int currentSharedPeak, Vector<SharedMS1PeakVO> sharedPeaks, Hashtable<Integer,Boolean> disentangledPeaks, SharedPeakContributionVO contr, Hashtable<String,LipidParameterSet> rtHits){
    Vector results = new Vector();
    Vector<LipidParameterSet> uniqueContributionsOfSpecies = areThereUniqueContributions(currentSharedPeak, sharedPeaks, disentangledPeaks, contr, rtHits);
    results.add(uniqueContributionsOfSpecies);
    if (uniqueContributionsOfSpecies.size()>0)
      results.add(getOtherSharedContributions(currentSharedPeak, sharedPeaks, disentangledPeaks, contr.getQuantVO()));
    return results;
  }

  
  /**
   * searches for unique contributions - peaks that belong only to one lipid species
   * @param currentSharedPeak the index of the currently shared peak under observation
   * @param sharedPeaks the vector containing all shared peak objects
   * @param disentangledPeaks boolean hash defining if shared peak instances were already disentangled and only one species is left
   * @param contr the current contributing partner where the unique and shared instances have to be found
   * @param rtHits all found hits currently valid hits
   * @return
   */
  private Vector<LipidParameterSet> areThereUniqueContributions(int currentSharedPeak, Vector<SharedMS1PeakVO> sharedPeaks, Hashtable<Integer,Boolean> disentangledPeaks, SharedPeakContributionVO contr, Hashtable<String,LipidParameterSet> rtHits){
    Vector<LipidParameterSet> uniqueOnes = new Vector<LipidParameterSet>();
    for (int i=0; i!=sharedPeaks.size();i++){
      if (i==currentSharedPeak) continue;
      if (!disentangledPeaks.get(i)) continue;
      //if it is an entangled peak, only this one should be left
      uniqueOnes.add(sharedPeaks.get(i).getPartners().get(0).getSet());
    }
    for (String rt : rtHits.keySet()){
      LipidParameterSet set = rtHits.get(rt);
      if (!(set instanceof LipidomicsMSnSet)) continue;
      if (!isHitASharedObject(contr.getQuantVO(), set, sharedPeaks)) uniqueOnes.add(set);
    }
    return uniqueOnes;
  }
  
  
  /**
   * checks if this LipidParameterSet is among the shared peak instances
   * @param anal the quantitation object holding information about the species
   * @param set the LipidParameterSet to be checked
   * @param sharedPeaks the vector containing all shared peak objects
   * @return
   */
  private boolean isHitASharedObject(QuantVO anal, LipidParameterSet set, Vector<SharedMS1PeakVO> sharedPeaks){
    for (SharedMS1PeakVO shared : sharedPeaks){
      for (SharedPeakContributionVO contr : shared.getPartners()){
        QuantVO other = contr.getQuantVO();
        if (anal.equals(other) && set.getRt().equalsIgnoreCase(contr.getSet().getRt())){
          return true;
        }
      }
    }
    return false;
  }
  
  
  /**
   * checks if there are other peaks from this species that are shared with other lipid species
   * @param currentSharedPeak the index of the currently shared peak under observation
   * @param sharedPeaks the vector containing all shared peak objects
   * @param disentangledPeaks boolean hash defining if shared peak instances were already disentangled and only one species is left
   * @param anal the quantitation object holding information about the species
   * @return a list of indices that name other shared instances that overlap as well
   */
  private Vector<Integer> getOtherSharedContributions (int currentSharedPeak, Vector<SharedMS1PeakVO> sharedPeaks, Hashtable<Integer,Boolean> disentangledPeaks, QuantVO anal){
    Vector<Integer> otherShared = new Vector<Integer>();
    for (int i=0; i!=sharedPeaks.size();i++){
      if (i==currentSharedPeak) continue;
      if (disentangledPeaks.containsKey(i) && disentangledPeaks.get(i)) continue;
      for (SharedPeakContributionVO contr : sharedPeaks.get(i).getPartners()){
        if (contr.getQuantVO().equals(anal)){
          otherShared.add(i);
          break;
        }
      }
    }
    return otherShared;
  }

  
  /**
   * detects if all of the unique detections are on one side of the contribution under observation
   * @param rt the retention time of the current contribution under observation
   * @param sets the other unique LipidParameterSets of this species
   * @return 0 if there is no direction; -1 if the peaks are in negative RT direction, +1 if the peaks are in positive RT direction
   */
  private int detectDirectionOfUniqueHits(double rt, Vector<LipidParameterSet> sets){
    int dir = 0;
    for (int i=0; i!=sets.size();i++){
      double other = Double.parseDouble(sets.get(i).getRt());
      if (i==0){
        if (other>rt) dir = 1;
        else dir = -1;
      }else{
        if (dir==1){
          if (other<rt){
            dir = 0;
            break;
          }
        }else{
          if (other>rt){
            dir = 0;
            break;
          }
        }
      }
    }
    return dir;
  }
  
  
  /**
   * returns all found indices of the shared objects that are higher than the currentIndex
   * @param currentIndex the index of the currently shared peak under observation
   * @param sharedInsts a list of indices that name other shared instances that overlap as well
   * @return all found indices of the shared objects that are higher than the currentIndex
   */
  private Vector<Integer> getSharedDescendingOrder(int currentIndex, Vector<Integer> sharedInsts){
    Vector<Integer> descending = new Vector<Integer>();
    for (Integer sharedIndex : sharedInsts){
      if (sharedIndex<=currentIndex) continue;
      int count = 0;
      for (int i=0; i!=descending.size(); i++){
        if (descending.get(i)<sharedIndex) break;
        count++;
      }
      if (count==descending.size()) descending.add(sharedIndex);
      else descending.add(count,sharedIndex);
    }
    return descending;
  }

  
  /**
   * checks if the current contribution under observation is sufficiently timely separated from the unique identifications
   * @param minRt the minimum RT time that is required to make this check
   * @param contr the current contributing partner to be checked
   * @param uniques the other unique LipidParameterSets of this species
   * @return true if the current contribution under observation is sufficiently timely separated from the unique identifications
   */
  private boolean isSufficentlyFarAway (float minRt, SharedPeakContributionVO contr, Vector<LipidParameterSet> uniques){
    boolean isOK = true;
    float rt = Float.parseFloat(contr.getSet().getRt());
    for (LipidParameterSet set : uniques){
      float otherRt = Float.parseFloat(set.getRt());
      if (Math.abs(rt-otherRt)<minRt){
        isOK = false;
        break;
      }
    }
    return isOK;
  }
 
  
  /**
   * calculates a percental split of the MS1 peak intensity according to the distinct MSn fragments of the sharing partners
   * @param shared the VO containing the shared MS1 peak and information about their sharing partners
   */
  private void calculatePercentualSplitValueAccordingToMSn(SharedMS1PeakVO shared){
    float sum = 0;
    Hashtable<QuantVO,Float> contrSums = new Hashtable<QuantVO,Float>();
    for (SharedPeakContributionVO contr : shared.getPartners()){
      if (contr.getOldSet()!=null) contr.setSet(contr.getOldSet());
      float contrSum = 0f;
      for (CgProbe fragment : contr.getDistinctFragments().values()) contrSum += fragment.Area;
      contrSums.put(contr.getQuantVO(), contrSum);
      sum += contrSum;
    }
    for (SharedPeakContributionVO contr : shared.getPartners()){
      contr.getSet().setPercentalSplit((contrSums.get(contr.getQuantVO())*100f)/sum);
    }
  }

  
  public Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> getPeaksBeforeSplit() {
    return this.peaksBeforeSplit_;
  }
}
