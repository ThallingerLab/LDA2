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

package at.tugraz.genome.lda;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import at.tugraz.genome.lda.exception.NoRuleException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.msn.MSnAnalyzer;
import at.tugraz.genome.lda.msn.RulesContainer;
import at.tugraz.genome.lda.msn.vos.MSnNamingVO;
import at.tugraz.genome.lda.msn.vos.SharedMS1PeakVO;
import at.tugraz.genome.lda.msn.vos.SharedPeakContributionVO;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.LipidomicsAnalyzer;
import at.tugraz.genome.lda.quantification.LipidomicsChromatogram;
import at.tugraz.genome.lda.quantification.LipidomicsDefines;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.utils.StringFloatVO;
import at.tugraz.genome.lda.vos.QuantVO;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.quantification.CgAreaStatus;
import at.tugraz.genome.maspectras.quantification.CgException;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.maspectras.utils.Calculator;
import at.tugraz.genome.voutils.GeneralComparator;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class SingleQuantThread extends Thread
{
  
  private LipidomicsAnalyzer analyzer_;
  private QuantVO quantSet_;
  private int msLevel_;
  /** should the MSn quantitation be performed before the MS1 quantitation */
  private boolean msnFirst_;
  private boolean finished_ = false;
  private String errorString_ = null;
  /** are there any MSn spectra present*/
  private boolean msnSpectraPresent_;
  private Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> results_;
  private Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> ms2RemovedHits_;
  /** if a peak split has to be removed because a split partner has a wrong retention time, the unsplit peak version is stored*/
  private Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> peaksBeforeSplit_;

  public SingleQuantThread(LipidomicsAnalyzer analyzer, QuantVO quantSet, int msLevel, boolean msnFirst){
    this.analyzer_ = analyzer;
    this.quantSet_ = quantSet;
    this.msLevel_ = msLevel;
    this.msnFirst_ = msnFirst;
    finished_ = false;
    errorString_ = null;
    msnSpectraPresent_ = false;
  }
  
  public void run(){
    try{
      ms2RemovedHits_ = new Hashtable<QuantVO,Hashtable<String,LipidParameterSet>>();
      peaksBeforeSplit_ = new Hashtable<QuantVO,Hashtable<String,LipidParameterSet>>();
      results_ = startSingleQuantification(analyzer_, quantSet_, msLevel_);
    } catch (Exception ex){
      ex.printStackTrace();
      errorString_ = ex.toString();
      
    }
    finished_ = true;
  }
  
  private Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> startSingleQuantification(LipidomicsAnalyzer analyzer, QuantVO quantSet, int msLevel) throws CgException {
    System.out.println("massOfInterest: "+quantSet.getAnalyteMass());
    Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>> isotopicProbes = null;
    float defaultRelativeAreaCutoff = analyzer.getRelativeAreaCutoff();
    float defaultRelativeFarAreaCutoff = analyzer.getRelativeFarAreaCutoff();
    int peakDiscardingAreaFactor = analyzer.getPeakDiscardingAreaFactor();
    Vector<QuantVO> quantVOs = new Vector<QuantVO>();
    quantVOs.add(quantSet);
    quantVOs.addAll(quantSet.getOtherIsobaricSpecies());
    String cutoff = null;
    for (QuantVO quant : quantVOs){
      ms2RemovedHits_.put(quant, new Hashtable<String,LipidParameterSet>());
      peaksBeforeSplit_.put(quant, new Hashtable<String,LipidParameterSet>());
      if (LipidomicsConstants.isMS2()){
        try{
          String aCutoff = RulesContainer.getMS1PeakCutoff(StaticUtils.getRuleName(quant.getAnalyteClass(),quant.getModName()));
          if (aCutoff!=null && Float.parseFloat(aCutoff)<analyzer.getRelativeFarAreaCutoff() && (cutoff==null || Float.parseFloat(aCutoff)<Float.parseFloat(cutoff))){
            cutoff = aCutoff;
          }
        } catch(NoRuleException nrx){
        } catch (RulesException e) {
          e.printStackTrace();
        } catch (IOException e) {
          e.printStackTrace();
        } catch (SpectrummillParserException e) {
          e.printStackTrace();
        }
      }
    }
    if (cutoff!=null){
      float cut = Float.parseFloat(cutoff);
      int peakDiscardingFactor = (int)(1f/cut);
      analyzer.setAreaCutoffs(cut, cut, peakDiscardingFactor, true);
    }

    
    if (msnFirst_){
      try {
        Vector<Float> rts = new Vector<Float>();
        MSnAnalyzer msnAnalyzer = null;
        for (QuantVO oneSet : quantVOs){
          msnAnalyzer = new MSnAnalyzer(oneSet.getAnalyteClass(),oneSet.getModName(),oneSet.getAnalyteMass(), (double)LipidomicsConstants.getMs2PrecursorTolerance(),
              oneSet.getAnalyteName(),oneSet.getDbs(),oneSet.getAnalyteFormula(), oneSet.getModFormula(),oneSet.getCharge(),analyzer,false,quantVOs.size()>1);
          rts.addAll(msnAnalyzer.getFoundMatchingSpectraTimes());
        }
        if (rts.size()>0){
          isotopicProbes = analyzer.processByMzProbabsAndPossibleRetentionTime((float)quantSet.getAnalyteMass(),quantSet.getCharge(),
              rts, quantSet.getUsedMinusTime(), quantSet.getUsedPlusTime(), LipidomicsDefines.MINUTES, 
              quantSet.getMustMatchProbabs(), quantSet.getProbabs(),msLevel,quantSet.getNegativeStartValue()<0);
        } else {
          if (msnAnalyzer!=null && msnAnalyzer.areMSnSpectraPresent()) msnSpectraPresent_ = true;
        }
      }
      catch (RulesException | IOException | SpectrummillParserException e) {
        e.printStackTrace();
      }

    } else {
      if (LipidomicsConstants.isShotgun()){
        isotopicProbes = analyzer.processShotgunData((float)quantSet.getAnalyteMass(),quantSet.getCharge(),msLevel,quantSet.getProbabs().size());
      } else if (quantSet.getIsobaricRetTime_()>0){
        if (LipidomicsConstants.use3D()){
          isotopicProbes = analyzer.processByMzProbabsAndPossibleRetentionTime((float)quantSet.getAnalyteMass(),quantSet.getCharge(),
            quantSet.getIsobaricRetTime_(), quantSet.getIsobaricMinusTime_(), quantSet.getUsedPlusTime(), LipidomicsDefines.MINUTES, 
            quantSet.getMustMatchProbabs(), quantSet.getProbabs(),msLevel,quantSet.getNegativeStartValue()<0);
        }else
          isotopicProbes = analyzer.processByMzAndRetentionTime((float)quantSet.getAnalyteMass(), quantSet.getCharge(),quantSet.getIsobaricRetTime_(), 
            quantSet.getIsobaricPlusTime_(), quantSet.getUsedPlusTime(), LipidomicsDefines.MINUTES,quantSet.getMustMatchProbabs(),
            quantSet.getProbabs(),msLevel,quantSet.getNegativeStartValue()<0);
      }else {
        if (LipidomicsConstants.use3D()){
          isotopicProbes = analyzer.processByMzProbabsAndPossibleRetentionTime((float)quantSet.getAnalyteMass(),quantSet.getCharge(),  
            -1f, -1f, -1f, -1, quantSet.getMustMatchProbabs(), quantSet.getProbabs(),msLevel,quantSet.getNegativeStartValue()<0);
        }else{
          isotopicProbes = analyzer.processByMzAndProbabs((float)quantSet.getAnalyteMass(),quantSet.getCharge(), quantSet.getMustMatchProbabs(),quantSet.getProbabs(),msLevel,quantSet.getNegativeStartValue()<0);
        }
      }
    }
    Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> hitsAccordingToQuant = new Hashtable<QuantVO,Hashtable<String,LipidParameterSet>>();
    //TODO: this is a quick hack - this can be improved by putting the MSnAanlyzer story in a separate method
    if (isotopicProbes!=null && isotopicProbes.size()>0){
      //first, generate the LipidParameterSet value objects
      Hashtable<Integer,Hashtable<QuantVO,LipidParameterSet>> isobarsOfAllSpecies = new Hashtable<Integer,Hashtable<QuantVO,LipidParameterSet>>();
      for (Integer key : isotopicProbes.keySet()){
        Hashtable<Integer,Vector<CgProbe>> oneHit = isotopicProbes.get(key);
        Hashtable<QuantVO,LipidParameterSet> isobars = new Hashtable<QuantVO,LipidParameterSet>();
        for (QuantVO oneSet : quantVOs){
          LipidParameterSet param = createLipidParameterSet(oneHit,oneSet.getNegativeStartValue(), (float)oneSet.getAnalyteMass(),
              oneSet.getAnalyteName(), oneSet.getDbs(), oneSet.getModName(), oneSet.getAnalyteFormula(), oneSet.getModFormula(), 
              oneSet.getCharge());
          if (LipidomicsConstants.isShotgun()){
            adaptMzValuesOfShotgunHits(param);
            isobars.put(oneSet, param);
          }else{
            String rt = param.getRt();
            float rtValue = Float.parseFloat(rt);
            if (oneSet.getRetTime()>0 && ((rtValue<(oneSet.getRetTime()-oneSet.getUsedMinusTime()))||(rtValue>(oneSet.getRetTime()+oneSet.getUsedPlusTime())))) continue;
            isobars.put(oneSet, param);
          }
        }
        if (isobars.size()>0) isobarsOfAllSpecies.put(key, isobars);
      }
      //second create m/z precursor and RT range VOs for each QuantVO for the caching of MS/MS spectra
      if (LipidomicsConstants.isMS2()){
        Hashtable<QuantVO,Hashtable<Integer,LipidParameterSet>> foundForQuantVO = new Hashtable<QuantVO,Hashtable<Integer,LipidParameterSet>>();
        for (QuantVO oneSet : quantVOs){
          float startMz = (float)oneSet.getAnalyteMass()-LipidomicsConstants.getMs2PrecursorTolerance();
          float stopMz = (float)oneSet.getAnalyteMass()+LipidomicsConstants.getMs2PrecursorTolerance();
          float lowestTime = Float.MAX_VALUE;
          float highestTime = 0f;
          if (LipidomicsConstants.isShotgun()){
            lowestTime = 0f;
            highestTime = Float.MAX_VALUE;
          }else{
            for (Hashtable<QuantVO,LipidParameterSet> isobars : isobarsOfAllSpecies.values()){
              if (!isobars.containsKey(oneSet)) continue;
              LipidParameterSet param = isobars.get(oneSet);
              float[] startStop = analyzer_.getStartStopTimeFromProbes(param.getIsotopicProbes().get(0));
              if (startStop[0]<lowestTime) lowestTime = startStop[0];
              if (startStop[1]>highestTime) highestTime = startStop[1];
            }
          }
          // third, prepare the spectral cache
          analyzer.prepareMSnSpectraCache(startMz, stopMz, lowestTime, highestTime);
          // fourth, do the MS2 detection
          Hashtable<Integer,LipidParameterSet> sameRt = new Hashtable<Integer,LipidParameterSet>();
          for (Integer key : isobarsOfAllSpecies.keySet()){
            Hashtable<QuantVO,LipidParameterSet> isobars = isobarsOfAllSpecies.get(key);
            if (!isobars.containsKey(oneSet)) continue;
            LipidParameterSet param = isobars.get(oneSet);
            boolean addHit = true;
            if (LipidomicsConstants.isMS2()){
              try {
                MSnAnalyzer msnAnalyzer = new MSnAnalyzer(oneSet.getAnalyteClass(),oneSet.getModName(),param,analyzer_,oneSet,false,quantVOs.size()>1);
                int msIdentOrder = RulesContainer.ORDER_MS1_FIRST;
                try{
                  msIdentOrder = RulesContainer.getMSIdentificationOrder(StaticUtils.getRuleName(oneSet.getAnalyteClass(),oneSet.getModName()));;
                } catch(Exception ex){
                }
                if (msnAnalyzer.checkStatus()==LipidomicsMSnSet.DISCARD_HIT || (msIdentOrder==RulesContainer.ORDER_MSN_ONLY && msnAnalyzer.checkStatus()<LipidomicsMSnSet.HEAD_GROUP_DETECTED)) addHit = false;
                else param = msnAnalyzer.getResult();
              }
              catch (RulesException e) {
                e.printStackTrace();
              }
              catch (IOException e) {
                e.printStackTrace();
              }
              catch (SpectrummillParserException e) {
                e.printStackTrace();
              }
            }
            if (addHit) sameRt.put(key, param);
            else {
              try{
                if (!LipidomicsConstants.isShotgun() && RulesContainer.isRtPostprocessing(StaticUtils.getRuleName(oneSet.getAnalyteClass(),oneSet.getModName())) && 
                    RulesContainer.correctRtForParallelModel(StaticUtils.getRuleName(oneSet.getAnalyteClass(),oneSet.getModName()))){
                  String rt = param.getRt();
                  ms2RemovedHits_.get(oneSet).put(rt, param);
                }
              } catch(NoRuleException nrx){
              } catch (RulesException e) {
                e.printStackTrace();
              } catch (IOException e) {
                e.printStackTrace();
              } catch (SpectrummillParserException e) {
                e.printStackTrace();
              }
            }
          }
          if (sameRt.size()>0)
            foundForQuantVO.put(oneSet, sameRt);
        }
        //fifth reorganize to group isobars together
        Hashtable<Integer,Hashtable<QuantVO,LipidParameterSet>> newIsobarsOfAllSpecies = new Hashtable<Integer,Hashtable<QuantVO,LipidParameterSet>>();
        for (QuantVO oneSet : foundForQuantVO.keySet()){
          Hashtable<Integer,LipidParameterSet> sameRt = foundForQuantVO.get(oneSet);
          for (Integer key : sameRt.keySet()){
            LipidParameterSet param = sameRt.get(key);
            Hashtable<QuantVO,LipidParameterSet> isobars = new Hashtable<QuantVO,LipidParameterSet>();
            if (newIsobarsOfAllSpecies.containsKey(key)) isobars = newIsobarsOfAllSpecies.get(key);
            isobars.put(oneSet, param);
            newIsobarsOfAllSpecies.put(key, isobars);
          }
        }
        isobarsOfAllSpecies = newIsobarsOfAllSpecies;
      }
      
      
      for (Hashtable<QuantVO,LipidParameterSet> isobarHitsOfOneSpecies : isobarsOfAllSpecies.values()){  
        // if there is only one found -> retry with absolute settings
        if (isobarHitsOfOneSpecies.size()==1){
          boolean addHit = true;
          QuantVO oneSet = isobarHitsOfOneSpecies.keySet().iterator().next();
          LipidParameterSet param = isobarHitsOfOneSpecies.get(oneSet);
          isobarHitsOfOneSpecies = new Hashtable<QuantVO,LipidParameterSet>();
//          System.out.println("!!!!!!! "+oneSet.getAnalyteClass()+param.getNameString()+"_"+oneSet.getModName()+" "+param.getRt());
          if (LipidomicsConstants.isMS2()){
            try {
              //TODO: the parameter before the last one is set to true in the meantime - maybe play around with caching in the future to improve calculation time
              MSnAnalyzer msnAnalyzer = new MSnAnalyzer(oneSet.getAnalyteClass(),oneSet.getModName(),param,analyzer_,oneSet,true,false);  
              if (msnAnalyzer.checkStatus()==LipidomicsMSnSet.DISCARD_HIT) addHit = false;
              else param = msnAnalyzer.getResult();
//              System.out.println("Status: "+msnAnalyzer.checkStatus());
            }
            catch (RulesException e) {
              e.printStackTrace();
            }
            catch (IOException e) {
              e.printStackTrace();
            }
            catch (SpectrummillParserException e) {
              e.printStackTrace();
            }
          }
          if (addHit) isobarHitsOfOneSpecies.put(oneSet,param);
          else{
            if (LipidomicsConstants.isMS2()){
              try{
                if (RulesContainer.isRtPostprocessing(StaticUtils.getRuleName(oneSet.getAnalyteClass(),oneSet.getModName())) && 
                    RulesContainer.correctRtForParallelModel(StaticUtils.getRuleName(oneSet.getAnalyteClass(),oneSet.getModName()))){
                  ms2RemovedHits_.get(oneSet).put(param.getRt(), param);
                }
              } catch(NoRuleException nrx){
              } catch (RulesException e) {
                e.printStackTrace();
              } catch (IOException e) {
                e.printStackTrace();
              } catch (SpectrummillParserException e) {
                e.printStackTrace();
              }
            }
          }
        }
        for (QuantVO quant : isobarHitsOfOneSpecies.keySet()){
          LipidParameterSet set = isobarHitsOfOneSpecies.get(quant);
          Hashtable<String,LipidParameterSet> hitsOfOneMod = new Hashtable<String,LipidParameterSet>();
          if (hitsAccordingToQuant.containsKey(quant)) hitsOfOneMod = hitsAccordingToQuant.get(quant);
          String rt = "";
          if (!LipidomicsConstants.isShotgun())
            rt = set.getRt();
          hitsOfOneMod.put(rt, set);
          hitsAccordingToQuant.put(quant, hitsOfOneMod);
        }
      }
    }
    //TODO: end of quick hack
    
    // this is for detection of peaks that are shared by different lipid classes
    if (quantVOs.size()>1){
      hitsAccordingToQuant = this.disentagleSharedMS1Peaks(hitsAccordingToQuant,analyzer,msLevel);
    }
    
    for (QuantVO oneSet : hitsAccordingToQuant.keySet()){
      Hashtable<String,LipidParameterSet> hitsOfOneMod = hitsAccordingToQuant.get(oneSet);
      if (LipidomicsConstants.isMS2() && hitsOfOneMod.size()>0){
        try {
          String unionTimeString = RulesContainer.getPeakUnionTime(StaticUtils.getRuleName(oneSet.getAnalyteClass(),oneSet.getModName()));
          if (unionTimeString != null && unionTimeString.length()>0){
            float unionTime = Float.parseFloat(unionTimeString);
            if (unionTime>0){
              boolean ignorePositionalEvidence = RulesContainer.isUnionWithoutPosition(StaticUtils.getRuleName(oneSet.getAnalyteClass(),oneSet.getModName()));
              hitsOfOneMod = uniteHits(analyzer,hitsOfOneMod,unionTime, oneSet,ignorePositionalEvidence);
              hitsAccordingToQuant.put(oneSet, hitsOfOneMod);
            }
          }
        }
        catch (RulesException | NoRuleException | IOException
            | SpectrummillParserException e) {
        }
      }      
    }
    
    if (LipidomicsConstants.isMS2()){
      if (cutoff!=null){
        analyzer.setAreaCutoffs(defaultRelativeAreaCutoff, defaultRelativeFarAreaCutoff, peakDiscardingAreaFactor, false);
      }
    }
    return hitsAccordingToQuant;
  }
  
  /**
   * this method splits the shared MS1 peaks according to the possibilities delivered by MS2 hits
   * @param hitsAccordingToQuant the results of the quantitation
   * @param analyzer the analyzer for quantifying hits
   * @param msLevel the level where the peak is oberlapping
   * @return the disentangled results
   * @throws CgException
   */
  @SuppressWarnings("unchecked")
  private Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> disentagleSharedMS1Peaks(Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> hitsAccordingToQuant, LipidomicsAnalyzer analyzer, int msLevel) throws CgException{
    //first detect which peaks are shared by different lipid classes
    Vector<SharedMS1PeakVO> sharedPeaks = detectSharedMS1PeakInstances(hitsAccordingToQuant);
    Vector<Integer> sharedToRemove = new Vector<Integer>();
//    System.out.println(sharedPeaks.size());
    for (int i=0; i!=sharedPeaks.size();i++){
      SharedMS1PeakVO shared = sharedPeaks.get(i);
      shared.checkForDistinctFragments();
      if (shared.hasAnyPartnerDistinctFragments()){
        removeDetectionsWhereDistinctFragmentsAreMissing(shared,hitsAccordingToQuant);
        if (shared.getPartners().size()==1){
          recheckMS2HitWithoutIsobar(shared,hitsAccordingToQuant);
          sharedToRemove.add(i);
        } else {
          int previousPartnerSize = shared.getPartners().size()+1;
          while (previousPartnerSize>shared.getPartners().size()){
            previousPartnerSize = shared.getPartners().size();
            shared.checkForDistinctFragments();
            removeDetectionsWhereSpectrumCoverageWithoutSharedIsToLow(shared,hitsAccordingToQuant,analyzer);
          }
          //System.out.println("!!!!!!!!!!!!!!!!!! "+shared.getPartners().size());
          
          if (shared.getPartners().size()==1){
            recheckMS2HitWithoutIsobar(shared,hitsAccordingToQuant);
            sharedToRemove.add(i);
          }
        }
      }
      if (!shared.hasAnyPartnerDistinctFragments()) {
        for (SharedPeakContributionVO contr : shared.getPartners()){
          contr.getSet().setPercentalSplit(100f);
        }
        sharedToRemove.add(i);
        if (shared.getPartners().size()>0) checkIfAllPartnerFragmentsFulfillSpectrumCoverage(shared,hitsAccordingToQuant,analyzer);
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
          removePartner(contr,shared,hitsAccordingToQuant);
          if (shared.getPartners().size()==1){
            recheckMS2HitWithoutIsobar(shared,hitsAccordingToQuant);
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
        if (!hitsAccordingToQuant.containsKey(contr.getQuantVO())) continue;
        @SuppressWarnings("rawtypes")
        Vector uniqueAndShared = getUniqueAndSharedContributions(i,sharedPeaks,disentangledPeaks,contr,hitsAccordingToQuant.get(contr.getQuantVO()));
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
              if (!hitsAccordingToQuant.containsKey(otherContr.getQuantVO())) continue;
              @SuppressWarnings("rawtypes")
              Vector otherUniqueAndShared = getUniqueAndSharedContributions(index,sharedPeaks,disentangledPeaks,otherContr,hitsAccordingToQuant.get(otherContr.getQuantVO()));
              if (otherUniqueAndShared.size()<2 || ((Vector<LipidParameterSet>)otherUniqueAndShared.get(0)).size()==0) continue;
              Vector<LipidParameterSet> otherUniques = (Vector<LipidParameterSet>)uniqueAndShared.get(0);
              int otherDirection = detectDirectionOfUniqueHits(Double.parseDouble(otherContr.getSet().getRt()), otherUniques);
              if (otherDirection!=direction) continue;
              if (isSufficentlyFarAway(rtDiff,otherContr,otherUniques) && other.isSpectrumContributionMuchLower(exclValue,otherContr)){
                removePartner(otherContr,other,hitsAccordingToQuant);
                if (other.getPartners().size()==1){
                  recheckMS2HitWithoutIsobar(other,hitsAccordingToQuant);
                  sharedToRemove.add(i);
                  disentangledPeaks.put(index, true);
                }
              }
            }
            
          }
          if (!hitsAccordingToQuant.containsKey(contr.getQuantVO())) continue;
          uniques = (Vector<LipidParameterSet>)getUniqueAndSharedContributions(i,sharedPeaks,disentangledPeaks,contr,hitsAccordingToQuant.get(contr.getQuantVO())).get(0);
        }
        if (isSufficentlyFarAway(rtDiff,contr,uniques) && shared.isSpectrumContributionMuchLower(exclValue,contr)){
          removePartner(contr,shared,hitsAccordingToQuant);
          if (shared.getPartners().size()==1){
            recheckMS2HitWithoutIsobar(shared,hitsAccordingToQuant);
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
      chrom = analyzer.readOneChromatogram((float)quant.getAnalyteMass(),msLevel);
    }
    for (SharedMS1PeakVO shared : sharedPeaks){
      float start = Float.MAX_VALUE;
      float stop = 0f;
      Vector<LipidParameterSet> sets = new Vector<LipidParameterSet>();
      for (SharedPeakContributionVO contr : shared.getPartners()) sets.add(contr.getSet());
      float[] borders = analyzer.calculatePercentualBorders(percentForBordersAbsolute,chrom,sets);
      start = borders[0];
      stop = borders[1];
      float[] bordersRelative = analyzer.calculatePercentualBorders(percentForBordersRelative,chrom,sets);
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
      if (shared.getPartners().size()==2 && analyzer.countMSnSpectraOfRegion(start, stop,2)>1){
        Vector<SharedPeakContributionVO> contrs = MSnAnalyzer.splitTwoIsobaricPeaks(analyzer,start,stop,startRelative,stopRelative,shared,msLevel);
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
            hitsAccordingToQuant.get(contrs.get(0).getQuantVO()).put(contrs.get(0).getSet().getRt(), contrs.get(0).getSet());
            peaksBeforeSplit_.get(contrs.get(0).getQuantVO()).put(contrs.get(0).getSet().getRt(), contrs.get(0).getOldSet());
          }
          if (contrs.get(1).getSet()!=null){
            hitsAccordingToQuant.get(contrs.get(1).getQuantVO()).put(contrs.get(1).getSet().getRt(), contrs.get(1).getSet());
            peaksBeforeSplit_.get(contrs.get(1).getQuantVO()).put(contrs.get(1).getSet().getRt(), contrs.get(1).getOldSet());
          }
        // no MS1 separation possible -> split according to MS2 intensities
        } else {
          calculatePercentualSplitValueAccordingToMSn(shared);
          for (SharedPeakContributionVO contr : shared.getPartners()){
            hitsAccordingToQuant.get(contr.getQuantVO()).put(contr.getSet().getRt(), contr.getSet());
          }

        }
      // here the peak intensity is split according to the distinct fragments TODO: has to be implemented
      } else {
        calculatePercentualSplitValueAccordingToMSn(shared);
//        System.out.println(classes+" "+rts);
        for (SharedPeakContributionVO contr : shared.getPartners()){
          hitsAccordingToQuant.get(contr.getQuantVO()).put(contr.getSet().getRt(), contr.getSet());
        }
      }
    }
    
    return hitsAccordingToQuant;
  }
  
  /**
   * this method unite hits if the MSn evidence if their evidence conforms, or at least does not contradict
   * @param analyzer the analyzer for calculating the united intensity
   * @param hits the hits which are checked for a union
   * @param unionTime the time, the hits for a union may be apart
   * @param quantSet the quantification settings VO
   * @param ignorePositionalEvidence ignore separations based on positional divergences
   * @return all of the hits, which include the potentially united results 
   * @throws CgException
   */
  @SuppressWarnings("unchecked")
  private Hashtable<String,LipidParameterSet> uniteHits(LipidomicsAnalyzer analyzer, Hashtable<String,LipidParameterSet> hits, float unionTime,
      QuantVO quantSet,boolean ignorePositionalEvidence) throws CgException{
    List<StringFloatVO> rts = new ArrayList<StringFloatVO>();
    for (String rt : hits.keySet())rts.add(new StringFloatVO(rt));
    Collections.sort(rts,new GeneralComparator("at.tugraz.genome.lda.utils.StringFloatVO", "getFloatValue", "java.lang.Float"));
    Vector<MSnNamingVO> candidates = new Vector<MSnNamingVO>();
    for (StringFloatVO vo : rts){
      LipidParameterSet ms1Set = hits.get(vo.getStringValue());
      if (!(ms1Set instanceof LipidomicsMSnSet)) continue;
      candidates.add(new MSnNamingVO((LipidomicsMSnSet)ms1Set,unionTime));
    }
    if (candidates.size()==0) return hits;
    Vector<MSnNamingVO> results = new Vector<MSnNamingVO>();
    int highestRegardedEvidence = LipidomicsMSnSet.POSITION_DETECTED;
    if (ignorePositionalEvidence) highestRegardedEvidence = LipidomicsMSnSet.FRAGMENTS_DETECTED;
    //for this procedure, the values must be sorted, otherwise, it won't work
    for (int status=highestRegardedEvidence; status>=LipidomicsMSnSet.HEAD_GROUP_DETECTED; status--){
      Vector<MSnNamingVO> toAnalyze = new Vector<MSnNamingVO>(results);
      Vector<MSnNamingVO> smallerStatus = new Vector<MSnNamingVO>();
      Set<Integer> combined = new HashSet<Integer>();
      if (results.size()==0) toAnalyze = new Vector<MSnNamingVO>(candidates);
      results = new Vector<MSnNamingVO>();
      for (int i=0; i!=toAnalyze.size(); i++){
        if (combined.contains(i)) continue;
        MSnNamingVO cand = toAnalyze.get(i);
        if (cand.getStatus()<status){
          smallerStatus.add(cand);
          continue;
        }
        //now I check if there is another candidate containing the same evidence
        for (int j=(i+1); j<toAnalyze.size(); j++){
          if (combined.contains(j)) continue;
          MSnNamingVO other = toAnalyze.get(j);
          if (other.getStatus()<status || !cand.insideRange(other) || !cand.hasSameEvidence(status,other)) continue;
          cand.mergeWithOther(other);
          combined.add(j);
        }
        // now I check if a combined RT covers any lower evidence candidates
        for (int j=(i+1); j<toAnalyze.size(); j++){
          if (combined.contains(j)) continue;
          MSnNamingVO other = toAnalyze.get(j);
          if (other.getStatus()<status || !cand.insideCombinedRtRange(other)) continue;
          cand.mergeWithOther(other);
          combined.add(j);          
        }      
        //now I check if there is a another candidate with contradicting Evidence
        for (int j=(i+1); j<toAnalyze.size(); j++){
          if (combined.contains(j)) continue;
          MSnNamingVO other = toAnalyze.get(j);
          if (other.getStatus()<status || !cand.insideRange(other)) continue;
          if (!cand.hasSameEvidence(status, other) && (status<LipidomicsMSnSet.POSITION_DETECTED || cand.hasSameEvidence(LipidomicsMSnSet.FRAGMENTS_DETECTED, other))){
            cand.setHigherLevelLimits(other);
          }
        }

        //now I check if this candidate is covered by an existing result - normally this should never be true
        boolean isCovered = false;
        for (MSnNamingVO result : results){
          if (result.insideCombinedRtRange(cand)){
            result.mergeWithOther(cand);
            System.out.println("ATTENTION: a result in the merging process is covered by an added one - this should not be the case!!!");
            break;
          }
        }
        if (isCovered) continue;
        //now I check if this new object covers existing results - normally this should never be true
        Vector<Integer> resultsToRemove = new Vector<Integer>();
        for (int j=0; j!=results.size(); j++){
          MSnNamingVO result = results.get(j);
          if (cand.insideCombinedRtRange(result)){
            cand.mergeWithOther(result);
            resultsToRemove.add(j);
            System.out.println("ATTENTION: a candidate covers an already added result - this should not be the case!!!");
            break;
          }
        }
        for (int j=resultsToRemove.size()-1; j!=-1; j--){
          results.remove(resultsToRemove.get(j));
        }
        results.add(cand);
      }
      results.addAll(smallerStatus);
    }
    
    // now combine the MS1 hits that fall inside the defined ranges
    Hashtable<String,LipidParameterSet> finalHits = new Hashtable<String,LipidParameterSet>();
    if (results.size()>0){
      Hashtable<MSnNamingVO,Hashtable<String,LipidParameterSet>> grouped = new  Hashtable<MSnNamingVO,Hashtable<String,LipidParameterSet>>();
      for (MSnNamingVO rangeVO : results) grouped.put(rangeVO, new Hashtable<String,LipidParameterSet>());
      for (String rtString : hits.keySet()){
        float rt = Float.parseFloat(rtString);
        LipidParameterSet set = hits.get(rtString);
        Vector<MSnNamingVO> insideRanges = new Vector<MSnNamingVO>();
        for (MSnNamingVO rangeVO : grouped.keySet()){
          if (rangeVO.insideRange(rt) && rangeVO.getMs1Name().equalsIgnoreCase(set.getNameStringWithoutRt())){
            insideRanges.add(rangeVO);
          }
        }
        if (insideRanges.size()==0) finalHits.put(rtString, set);
        else {
          MSnNamingVO rangeVO = detectClosestRange(Float.parseFloat(rtString), insideRanges);
          Hashtable<String,LipidParameterSet> toMerge = grouped.get(rangeVO);
          toMerge.put(rtString, set);
          grouped.put(rangeVO, toMerge);
        }
      }
      for (Hashtable<String,LipidParameterSet> forMerging : grouped.values()){
        if (forMerging.size()==1){
          String rt = forMerging.keySet().iterator().next();
          finalHits.put(rt, forMerging.get(rt));
          continue;
        }
        LipidParameterSet set = createLipidParameterSet(analyzer.mergeIdentifications(forMerging.values()),quantSet.getNegativeStartValue(),
            (float)quantSet.getAnalyteMass(), quantSet.getAnalyteName(), quantSet.getDbs(), quantSet.getModName(), quantSet.getAnalyteFormula(),
            quantSet.getModFormula(), quantSet.getCharge());
        try {
          //TODO: the parameter before the last one is set to true in the meantime - maybe play around with caching in the future to improve calculation time
          MSnAnalyzer msnAnalyzer = new MSnAnalyzer(quantSet.getAnalyteClass(),quantSet.getModName(),set,analyzer_,quantSet,true,false);  
          if (msnAnalyzer.checkStatus()!=LipidomicsMSnSet.DISCARD_HIT) set = msnAnalyzer.getResult();
        }
        catch (RulesException e) {
          e.printStackTrace();
        }
        catch (IOException e) {
          e.printStackTrace();
        }
        catch (SpectrummillParserException e) {
          e.printStackTrace();
        }
        finalHits.put(set.getRt(), set);
        
      }
    } else finalHits = new Hashtable<String,LipidParameterSet>(hits);
    return finalHits;
  }
  
  /**
   * method for creating a LipidParameterSet out of the detected peaks from the analyzer
   * @param oneHit the peaks of one identification
   * @param negativeStartValue is the isotopic series negative?
   * @param analyteMass the theoretical mass of the identification
   * @param analyteName the name of the analyte
   * @param dbs the number of the double bonds of the analyte
   * @param modName the name of the adduct/deduct
   * @param analyteFormula the formula of the analyte
   * @param modFormula the name of the modification of the adduct/deduct
   * @param charge the charge state of the identification
   * @return LipidParameterSet object 
   */
  private LipidParameterSet createLipidParameterSet(Hashtable<Integer,Vector<CgProbe>> oneHit, int negativeStartValue,
      float analyteMass, String analyteName, int dbs, String modName, String analyteFormula, String modFormula, int charge){
    Vector<Vector<CgProbe>> isotopicProbes2 = new Vector<Vector<CgProbe>>();
    // k is the isotope number
    float totalArea = 0;
    String rt = "";
    for (int k=0;k!=oneHit.size();k++){
      Vector<CgProbe> probes = oneHit.get(k);
      if (probes !=null){
        isotopicProbes2.add(probes);
        Vector<Double> rts = new Vector<Double>();
        for (int j=0; j!=probes.size(); j++){
          CgProbe probe = probes.get(j);
          if (negativeStartValue<0)
            probe.isotopeNumber *= -1;
          if (probe!=null&&probe.AreaStatus==CgAreaStatus.OK){
            totalArea += probe.Area;
            if (!LipidomicsConstants.isShotgun())
              rts.add((double)probe.Peak);
          }
        }
        if (!LipidomicsConstants.isShotgun() && k==0){
          rt = Calculator.FormatNumberToString(Calculator.mean(rts)/60d,2d);
        }
      }
    }
    if (LipidomicsConstants.isShotgun())
      rt = null;
    LipidParameterSet param = new LipidParameterSet(analyteMass, analyteName, dbs, modName, rt, analyteFormula, modFormula,charge);
    param.LowerMzBand = LipidomicsConstants.getCoarseChromMzTolerance(analyteMass);
    param.UpperMzBand = LipidomicsConstants.getCoarseChromMzTolerance(analyteMass);
    param.Area = totalArea;
    param.setIsotopicProbes(isotopicProbes2);
    return param;
  }
  
  
  public boolean finished(){
    return this.finished_;
  }

  public String getErrorString(){
    return this.errorString_;
  }
  
  public Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> getResults(){
    return this.results_;
  }

  public Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> getMs2Removedhits(){
    return this.ms2RemovedHits_;
  }

  /**
   * 
   * @return whether there are any MSn spectra found for this precursor ion m/z value
   */
  public boolean areMSnSpectraPresent()
  {
    return msnSpectraPresent_;
  }

  /**
   * 
   * @return if this quantitation was based on a preceding MSn search
   */
  public boolean isMsnFirst()
  {
    return msnFirst_;
  }

  /**
   * 
   * @return the value object containing the quantitation parameters
   */
  public QuantVO getQuantSet()
  {
    return quantSet_;
  }
  
  
  /**
   * detects the closest rangeVO out of a set of given rangeVOs
   * @param rt the retention time of the hit, where the closest range shall be detected
   * @param ranges the rangeVOs
   * @return the closest RangeVO
   */
  private MSnNamingVO detectClosestRange(float rt, Vector<MSnNamingVO> ranges){
    float closest = Float.MAX_VALUE;
    MSnNamingVO close = null;
    for (MSnNamingVO range : ranges){
      float dist = range.absoluteDistanceToHardRange(rt);
      if (dist<closest){
        closest = dist;
        close = range;
      }
    }
    return close;
  }

  
  /**
   * 
   * @param hitsAccordingToQuant quantitation hits that are assigned to QuantVOs (first key); second key is the retention time
   * @return Vector of SharedMS1PeakVO which are peaks that match to more than one analyte, and contain SharedPeakContributionVO
   *         which hold the QuantVO, the LipidParameterSet, and information about distinct fragments
   */
  private Vector<SharedMS1PeakVO> detectSharedMS1PeakInstances(Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> hitsAccordingToQuant){
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
        sets.remove(partner.getSet().getRt());
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
        sets.remove(partner.getSet().getRt());
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
   * @throws CgException exception that is thrown when there is something wrong with the fragment detection
   */
  private void recheckMS2HitWithoutIsobar(SharedMS1PeakVO shared, Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> hitsAccordingToQuant) throws CgException{
    for (SharedPeakContributionVO contr : shared.getPartners()){
      QuantVO oneSet = contr.getQuantVO();
      Hashtable<String,LipidParameterSet> quantsOfMod = hitsAccordingToQuant.get(oneSet);
      try {
        //TODO: this is set to true in the meantime - maybe play around with caching in the future to improve calculation time
        MSnAnalyzer msnAnalyzer = new MSnAnalyzer(oneSet.getAnalyteClass(),oneSet.getModName(),contr.getSet(),analyzer_,oneSet,true,false);  
        if (msnAnalyzer.checkStatus()==LipidomicsMSnSet.DISCARD_HIT) quantsOfMod.remove(contr.getSet().getRt());
        else quantsOfMod.put(contr.getSet().getRt(),msnAnalyzer.getResult());
      }
      catch (RulesException e) {
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
        hitsAccordingToQuant.get(contr.getQuantVO()).remove(contr.getSet().getRt());
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
          hitsAccordingToQuant.get(cont.getQuantVO()).remove(cont.getSet().getRt());
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

  /**
   * if a peak split has to be removed because a split partner has a wrong retention time, the unsplit peak version is stored
   * @return the unsplit peak versions
   */
  public Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> getPeaksBeforeSplit()
  {
    return peaksBeforeSplit_;
  }
  
  /**
   * sets corresponding m/z values for the shotgun hits
   * @param param the LipidParameterSet
   */
  private void adaptMzValuesOfShotgunHits(LipidParameterSet param){
    Vector<Vector<CgProbe>> isoProbes = param.getIsotopicProbes();
    Vector<Vector<CgProbe>> correctedProbes = new Vector<Vector<CgProbe>>();
    float mzTolerance = 0f;
    for (int i=0; i!=isoProbes.size(); i++){
      Vector<CgProbe> probes = isoProbes.get(i);
      Vector<CgProbe> corrected = new Vector<CgProbe>();
      float mz = param.Mz[0]+i*LipidomicsConstants.getNeutronMass()/(float)param.getCharge();
      for (CgProbe aProbe : probes){
        mzTolerance = aProbe.LowerMzBand;
        CgProbe probe = new CgProbe(aProbe);
        probe.Mz = mz;
        corrected.add(probe);
      }
      correctedProbes.add(corrected);
    }
    param.setIsotopicProbes(correctedProbes);
    param.LowerMzBand = mzTolerance;
    param.UpperMzBand = mzTolerance;
  }
  
  
}
