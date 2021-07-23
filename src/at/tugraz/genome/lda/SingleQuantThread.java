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

import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.HydroxylationEncodingException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.exception.NoRuleException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.msn.MSnAnalyzer;
import at.tugraz.genome.lda.msn.MSnPeakSeparator;
import at.tugraz.genome.lda.msn.RulesContainer;
import at.tugraz.genome.lda.msn.vos.MSnNamingVO;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.LipidomicsAnalyzer;
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
  
  private Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> startSingleQuantification(LipidomicsAnalyzer analyzer, QuantVO quantSet, int msLevel) throws CgException, LipidCombinameEncodingException {
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
          msnAnalyzer = new MSnAnalyzer(oneSet.getAnalyteClass(),oneSet.getModName(),oneSet.getAnalyteMass(),(double)LipidomicsConstants.getMs2PrecursorTolerance(),
              oneSet.getAnalyteName(),oneSet.getDbs(),oneSet.getOhNumber(),oneSet.getAnalyteFormula(), oneSet.getModFormula(),oneSet.getCharge(),analyzer,false,
              quantVOs.size()>1);
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
      catch (RulesException | IOException | SpectrummillParserException | HydroxylationEncodingException | ChemicalFormulaException e) {
        e.printStackTrace();
      }

    } else {
      if (LipidomicsConstants.isShotgun()==LipidomicsConstants.SHOTGUN_TRUE){
        isotopicProbes = analyzer.processShotgunData((float)quantSet.getAnalyteMass(),quantSet.getCharge(),msLevel,quantSet.getProbabs().size());
      } else if (LipidomicsConstants.isShotgun()==LipidomicsConstants.SHOTGUN_PRM){
        //TODO: the MS-level is here always set to 2, and the ohNumber is always set to 0
        isotopicProbes = analyzer.processPrmData((float)quantSet.getAnalyteMass(),quantSet.getCharge(),2,quantSet.getAnalyteClass(),
            quantSet.getModName(), StaticUtils.generateLipidNameString(quantSet.getAnalyteName(), quantSet.getDbs(), -1),quantSet.getAnalyteFormula(),0);
//        for (Integer hitNumber : isotopicProbes.keySet()){
//          for (Integer isoNr: isotopicProbes.get(hitNumber).keySet()){
//            for (CgProbe probe : isotopicProbes.get(hitNumber).get(isoNr)){
//              System.out.println("! "+probe.Peak/60f+" "+probe.Area);
//            }
//          }
//        }
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
              oneSet.getAnalyteName(), oneSet.getDbs(), oneSet.getOhNumber(), oneSet.getModName(), oneSet.getAnalyteFormula(), oneSet.getModFormula(), 
              oneSet.getCharge());
          if (LipidomicsConstants.isShotgun()==LipidomicsConstants.SHOTGUN_TRUE){
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
          if (LipidomicsConstants.isShotgun()==LipidomicsConstants.SHOTGUN_TRUE){
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
              catch (HydroxylationEncodingException | ChemicalFormulaException | LipidCombinameEncodingException e) {
                e.printStackTrace();
              }
            }
            if (addHit) sameRt.put(key, param);
            else {
              try{
                if (LipidomicsConstants.isShotgun()!=LipidomicsConstants.SHOTGUN_TRUE && RulesContainer.isRtPostprocessing(StaticUtils.getRuleName(oneSet.getAnalyteClass(),oneSet.getModName())) && 
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
          if (addHit) isobarHitsOfOneSpecies.put(oneSet,param);
          else{
            if (LipidomicsConstants.isMS2()){
              try{
                if (LipidomicsConstants.isShotgun()!=LipidomicsConstants.SHOTGUN_TRUE && RulesContainer.isRtPostprocessing(StaticUtils.getRuleName(oneSet.getAnalyteClass(),oneSet.getModName())) && 
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
          if (LipidomicsConstants.isShotgun()!=LipidomicsConstants.SHOTGUN_TRUE)
            rt = set.getRt();
          hitsOfOneMod.put(rt, set);
          hitsAccordingToQuant.put(quant, hitsOfOneMod);
        }
      }
    }
    //TODO: end of quick hack
    
    // this is for detection of peaks that are shared by different lipid classes
    if (quantVOs.size()>1){
      Set<String> adductsThatRequireOtherAdducts = new HashSet<String>();
      String rule;
      for (QuantVO quant : hitsAccordingToQuant.keySet()) {
        try{
          rule = StaticUtils.getRuleName(quant.getAnalyteClass(),quant.getModName());
          if (RulesContainer.requiresOtherValidAdduct(rule))
            adductsThatRequireOtherAdducts.add(rule);
        } catch(NoRuleException nrx){
        } catch (RulesException e) {
          e.printStackTrace();
        } catch (IOException e) {
          e.printStackTrace();
        } catch (SpectrummillParserException e) {
          e.printStackTrace();
        }
      }
      //System.out.println("areHitsPresentThatRelyOnOtherAdducts: "+areHitsPresentThatRelyOnOtherAdducts);
      //this has to be done for hits where there are other adducts need to be detected - they have to be entagled later on
//      if (!areHitsPresentThatRelyOnOtherAdducts) {
      MSnPeakSeparator separator = new MSnPeakSeparator(hitsAccordingToQuant, peaksBeforeSplit_, analyzer, msLevel, adductsThatRequireOtherAdducts);
      hitsAccordingToQuant = separator.disentagleSharedMS1Peaks();
//      }
    }
    
    for (QuantVO oneSet : hitsAccordingToQuant.keySet()){
      Hashtable<String,LipidParameterSet> hitsOfOneMod = hitsAccordingToQuant.get(oneSet);
      if (LipidomicsConstants.isShotgun()!=LipidomicsConstants.SHOTGUN_TRUE && LipidomicsConstants.isMS2() && hitsOfOneMod.size()>0){
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
            | SpectrummillParserException | LipidCombinameEncodingException e) {
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
   * this method unite hits if the MSn evidence if their evidence conforms, or at least does not contradict
   * @param analyzer the analyzer for calculating the united intensity
   * @param hits the hits which are checked for a union
   * @param unionTime the time, the hits for a union may be apart
   * @param quantSet the quantification settings VO
   * @param ignorePositionalEvidence ignore separations based on positional divergences
   * @return all of the hits, which include the potentially united results 
   * @throws CgException
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  @SuppressWarnings({ "unchecked", "unlikely-arg-type" })
  private Hashtable<String,LipidParameterSet> uniteHits(LipidomicsAnalyzer analyzer, Hashtable<String,LipidParameterSet> hits, float unionTime,
      QuantVO quantSet,boolean ignorePositionalEvidence) throws CgException, LipidCombinameEncodingException{
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
            (float)quantSet.getAnalyteMass(), quantSet.getAnalyteName(), quantSet.getDbs(), quantSet.getOhNumber(), quantSet.getModName(),
            quantSet.getAnalyteFormula(), quantSet.getModFormula(), quantSet.getCharge());
        try {
          //TODO: the parameter before the last one is set to true in the meantime - maybe play around with caching in the future to improve calculation time
          MSnAnalyzer msnAnalyzer = new MSnAnalyzer(quantSet.getAnalyteClass(),quantSet.getModName(),set,analyzer_,quantSet,true,false);  
          if (msnAnalyzer.checkStatus()!=LipidomicsMSnSet.DISCARD_HIT) set = msnAnalyzer.getResult();
        }
        catch (RulesException | HydroxylationEncodingException | ChemicalFormulaException e) {
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
   * @param ohNumber the number of hydroxylation sites
   * @param modName the name of the adduct/deduct
   * @param analyteFormula the formula of the analyte
   * @param modFormula the name of the modification of the adduct/deduct
   * @param charge the charge state of the identification
   * @return LipidParameterSet object 
   */
  private LipidParameterSet createLipidParameterSet(Hashtable<Integer,Vector<CgProbe>> oneHit, int negativeStartValue,
      float analyteMass, String analyteName, int dbs, int ohNumber, String modName, String analyteFormula, String modFormula,
      int charge){
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
            if (LipidomicsConstants.isShotgun()!=LipidomicsConstants.SHOTGUN_TRUE)
              rts.add((double)probe.Peak);
          }
        }
        if (LipidomicsConstants.isShotgun()!=LipidomicsConstants.SHOTGUN_TRUE && k==0){
          rt = Calculator.FormatNumberToString(Calculator.mean(rts)/60d,2d);
        }
      }
    }
    if (LipidomicsConstants.isShotgun()==LipidomicsConstants.SHOTGUN_TRUE)
      rt = null;
    LipidParameterSet param = new LipidParameterSet(analyteMass, analyteName, dbs, ohNumber, modName, rt, analyteFormula, modFormula,charge);
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

  
  /**
   * @return the LipidomicsAnalyzer that was used for data access
   */
  public LipidomicsAnalyzer getAnalyzer()
  {
    return analyzer_;
  }
  
  
  
  
}
