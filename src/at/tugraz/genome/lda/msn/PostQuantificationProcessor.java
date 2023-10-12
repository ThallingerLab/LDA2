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

package at.tugraz.genome.lda.msn;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.exception.LMException;
import at.tugraz.genome.lda.exception.NoRuleException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.msn.parser.FragRuleParser;
import at.tugraz.genome.lda.msn.vos.RtPredictVO;
import at.tugraz.genome.lda.msn.vos.SharedMS1PeakVO;
import at.tugraz.genome.lda.msn.vos.SharedPeakContributionVO;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.utils.LMAsymptDecayThreeVariables;
import at.tugraz.genome.lda.utils.LMAsymptDecayTwoVariables;
import at.tugraz.genome.lda.utils.LMAsymptVarDecayThreeVariables;
import at.tugraz.genome.lda.utils.LMAsymptVarDecayTwoVariables;
import at.tugraz.genome.lda.utils.LMLogDecayThreeVariables;
import at.tugraz.genome.lda.utils.LMLogDecayTwoVariables;
import at.tugraz.genome.lda.utils.RangeInteger;
import at.tugraz.genome.lda.utils.LevenbergMarquardtOptimizer;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.QuantVO;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.maspectras.utils.Calculator;
import at.tugraz.genome.util.FloatMatrix;

/**
 * This class makes post processing operation after MS identification and quantification
 * @author Juergen Hartler
 *
 */
public class PostQuantificationProcessor
{
  
  /** the analytes that pass the filter*/
  private Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> results_;
  /** analytes that were removed before (usually by MSn) and serve as negative filter*/
  private Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> ms2Removed_;
  /** should the RT filter be calculated based on all modifications */
  private Hashtable<String,Boolean> adductInsensitiveRtFilter_;
  
  private final static int MINIMUM_MEASUREMENTS = 7;
  private final static float COUNTER_SERIES_TOLERANCE = 3f;
  /** the predicted models - key class name or rule name, depending on adduct insensitive filter*/
  private Hashtable<String,LevenbergMarquardtOptimizer> predictedModels_;
  /** does the predicted model respect OH*/
  private Hashtable<String,Boolean> respectOhs_;
  
  
  /**
   * constructor requiring the analytes that have to be filtered and the ones for the negative filter (the ones for the negative filter may be empty)
   * @param results analytes that have to pass the filter
   * @param ms2Removed analytes for the negative filter (may be empty)
   */
  public PostQuantificationProcessor(Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> results,
      Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> ms2Removed, Hashtable<String,Boolean> adductInsensitiveRtFilter){
    results_ = results;
    ms2Removed_ = ms2Removed;
    adductInsensitiveRtFilter_ = adductInsensitiveRtFilter;
    predictedModels_ = new Hashtable<String,LevenbergMarquardtOptimizer>();
    respectOhs_ = new Hashtable<String,Boolean>();
  }
  
  /**
   * starts the filtering process and returns the filtered data
   * @return filtered data
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the library is not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  public Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> processData() throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    results_ = correctByRetentionTimeSeries(results_,ms2Removed_,adductInsensitiveRtFilter_);
    return results_;
  }
  
  /**
   * predicts the RT of lipids where no MSn spectra are present - for calculation at these retnetion times
   * @param ms1ToPredict lipids to quantify where no MSn spectra are present
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the library is not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  @SuppressWarnings("unchecked")
  public Hashtable<String,Hashtable<String,RtPredictVO>> predictRetentionTimesBasedOnResults(Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> ms1ToPredict, Hashtable<String,Hashtable<String,RtPredictVO>> prevPredictions) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    Vector<Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>>> dataOrdered = orderData(results_,ms2Removed_,null);
    Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> postProcessData = dataOrdered.get(1);
    Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> negativeExamples = dataOrdered.get(2);
    Hashtable<String,Hashtable<String,RtPredictVO>> rtPredictVOs = new Hashtable<String,Hashtable<String,RtPredictVO>>();
    
    for (String className : ms1ToPredict.keySet()){
      if (!postProcessData.containsKey(className)) continue;
      Hashtable<String,RtPredictVO> classRtPredictVOs = new Hashtable<String,RtPredictVO>();
      for (String mod : postProcessData.get(className).keySet()){
        Hashtable<String,Hashtable<String,LipidParameterSet>> negatives = new Hashtable<String,Hashtable<String,LipidParameterSet>>();
        if (negativeExamples.containsKey(className) && negativeExamples.get(className).containsKey(mod)) negatives = negativeExamples.get(className).get(mod);
        String ruleName = StaticUtils.getRuleName(className, mod);
        Hashtable<String,Hashtable<String,LipidParameterSet>> result = postProcessData.get(className).get(mod);
        try {
          Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> paramsOrdered = new Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>();
          Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> negativesOrdered = new Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>();
          @SuppressWarnings("rawtypes")
          Vector ranges;
          RangeInteger cAtomsRange = null;
          Hashtable<Integer,RangeInteger> dbsRanges = null;
//          if (prevPredictions!=null && prevPredictions.containsKey(className) && prevPredictions.get(className).containsKey(mod)){
//            RtPredictVO predVO = prevPredictions.get(className).get(mod);
//            cAtomsRange = predVO.getcAtomsRange();
//            dbsRanges = predVO.getDbsRanges();
//          }else{
            ranges = groupAccordingToCAtomsAndDoubleBonds(result,new Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>(),paramsOrdered,className,ruleName);
            cAtomsRange = (RangeInteger)ranges.get(0);
            dbsRanges = (Hashtable<Integer,RangeInteger>)ranges.get(1);
//          }
          @SuppressWarnings({ "rawtypes", "unused" })
          Vector rangesNeg = groupAccordingToCAtomsAndDoubleBonds(negatives,new Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>(),negativesOrdered,className,ruleName);
          float tolerance = 4f;
          float maxDev = -1f;
          if (RulesContainer.getRetentionTimeMaxDeviation(ruleName)!=null) maxDev = new Float(RulesContainer.getRetentionTimeMaxDeviation(ruleName));
          
          //this is for oh
          boolean diffOh = false;
          int oh = -1;
          for (Hashtable<Integer,Hashtable<String,LipidParameterSet>> sameC : paramsOrdered.values()){
            for (Hashtable<String,LipidParameterSet> sameDbs : sameC.values()) {
              for (LipidParameterSet set : sameDbs.values()) {
                if (oh==-1)
                  oh = set.getOhNumber();
                else if (oh!=set.getOhNumber())
                  diffOh = true;
              }
            }
          }
          @SuppressWarnings("rawtypes")
          Vector resultsLM = doIterativeLMOptimization(paramsOrdered,paramsOrdered,cAtomsRange,dbsRanges,negativesOrdered,tolerance,maxDev,diffOh);
          LevenbergMarquardtOptimizer optimizer = (LevenbergMarquardtOptimizer)resultsLM.get(0);
          LevenbergMarquardtOptimizer counterModel = (LevenbergMarquardtOptimizer)resultsLM.get(1);
                    
          // the model was fitted - now check for which hits shall we make the prediction
          Hashtable<String,Hashtable<String,LipidParameterSet>> unprocResult = new Hashtable<String,Hashtable<String,LipidParameterSet>>(result); 
          for (String analyteName : ms1ToPredict.get(className).keySet()){
            Hashtable<String,QuantVO> analytesMod = ms1ToPredict.get(className).get(analyteName);
            if (!analytesMod.containsKey(mod) || unprocResult.containsKey(analyteName)) continue;
            QuantVO quantVO = analytesMod.get(mod);
            LipidParameterSet setForPred = new LipidParameterSet((float)quantVO.getAnalyteMass(), quantVO.getAnalyteName(),
                quantVO.getDbs(), quantVO.getModName(), -1.0, quantVO.getAnalyteFormula(), quantVO.getModFormula(),
                quantVO.getCharge(), quantVO.getOhNumber());
            Hashtable<String,LipidParameterSet> forPred = new Hashtable<String,LipidParameterSet>();
            forPred.put(setForPred.getRt(), setForPred);
            unprocResult.put(analyteName, forPred);
          }
          Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> unprocessed = new Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>();
          ranges = groupAccordingToCAtomsAndDoubleBonds(unprocResult,new Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>(),unprocessed,className,ruleName);
          RangeInteger cAtomsMaxRange = (RangeInteger)ranges.get(0);
          Hashtable<Integer,RangeInteger> dbsMaxRanges = (Hashtable<Integer,RangeInteger>)ranges.get(1);
          if (prevPredictions!=null && prevPredictions.containsKey(className) && prevPredictions.get(className).containsKey(mod)){
            RtPredictVO predVO = prevPredictions.get(className).get(mod);
            cAtomsRange = predVO.getcAtomsRange();
            dbsRanges = predVO.getDbsRanges();
          }
          @SuppressWarnings("rawtypes")
          Vector proposedRanges = proposeLMFilterRanges(unprocessed, cAtomsMaxRange, dbsMaxRanges, cAtomsRange, dbsRanges);
          RangeInteger newCAtomsRange = (RangeInteger)proposedRanges.get(0);
          Hashtable<Integer,RangeInteger> newDbsRanges = (Hashtable<Integer,RangeInteger>) proposedRanges.get(1);
//        this is only for debug purposes
//          for (int i=newCAtomsRange.getStart(); i<=newCAtomsRange.getStop(); i++){
//            if (!newDbsRanges.containsKey(i)) continue;
//            RangeInteger dbsRange = newDbsRanges.get(i);
//            System.out.println("Proposed Range calculation: "+i+":"+dbsRange.getStart()+"-"+dbsRange.getStop());
//          }
          RtPredictVO predVO = new RtPredictVO(optimizer,counterModel,newCAtomsRange,newDbsRanges);
          classRtPredictVOs.put(mod, predVO);
          
          Pattern cAtomsPattern =  Pattern.compile(RulesContainer.getCAtomsFromNamePattern(ruleName));
          Pattern dbsPattern =  Pattern.compile(RulesContainer.getDoubleBondsFromNamePattern(ruleName));
          for (String analyteName : ms1ToPredict.get(className).keySet()){
            Hashtable<String,QuantVO> analytesMod = ms1ToPredict.get(className).get(analyteName);
            if (!analytesMod.containsKey(mod)) continue;
            QuantVO quantVO = analytesMod.get(mod);
            Matcher cAtomsMatcher = cAtomsPattern.matcher(analyteName.split(LipidomicsConstants.CHAIN_MOD_SEPARATOR)[0]);
            if (!cAtomsMatcher.matches()) throw new RulesException("The analyte "+analyteName+" does not match the "+FragRuleParser.GENERAL_CATOMS_PARSE+" pattern \""+RulesContainer.getCAtomsFromNamePattern(ruleName)+"\" of the class "+ruleName+"!");
            Matcher dbsMatcher = dbsPattern.matcher(analyteName.split(LipidomicsConstants.CHAIN_MOD_SEPARATOR)[0]);
            int cAtoms = Integer.parseInt(cAtomsMatcher.group(1));
            if (!dbsMatcher.matches()) throw new RulesException("The analyte "+analyteName+" does not match the "+FragRuleParser.GENERAL_DBOND_PARSE+" pattern \""+RulesContainer.getDoubleBondsFromNamePattern(ruleName)+"\" of the class "+ruleName+"!");
            int dbs = Integer.parseInt(dbsMatcher.group(1));
// the lines with the 4/ are necessary if the prediction is based on consecutive model predictions
////            if (!newDbsRanges.containsKey(cAtoms)) continue;
////            RangeInteger dbsRange = newDbsRanges.get(cAtoms);
////            if (dbsRange.insideRange(dbs))
              quantVO.setRetTime(optimizer.calculateFitValue(new float[]{cAtoms,dbs}));
          }          
        }
        catch (LMException e) {
          System.out.println("Warning: "+className+"_"+mod+" was not RT filtered: "+e.getMessage());
        }
      }
      if (classRtPredictVOs.size()>0) rtPredictVOs.put(className, classRtPredictVOs);
    }
    return rtPredictVOs;
  }

  
  /**
   * the filter fits a curve to the hits that were identified by MSn and uses negative hits for a counter curve
   * this method groups the data according to lipid class and modification
   * @param unprocessed the analytes that need to be processed
   * @param ms2Removed hits removed by MSn - useable for counter curve
   * @return filtered the filtered analytes
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the library is not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   */
  private Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> correctByRetentionTimeSeries(Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> unprocessed,
      Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> ms2Removed, Hashtable<String,Boolean> adductInsensitiveRtFilter_) throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    Vector<Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>>> dataOrdered = orderData(unprocessed,ms2Removed,adductInsensitiveRtFilter_);
    Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> results = dataOrdered.get(0);
    Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> postProcessData = dataOrdered.get(1);
    Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> negativeExamples = dataOrdered.get(2);
    // start post processing for every class mod combination
    for (String className : postProcessData.keySet()){
      if (adductInsensitiveRtFilter_.get(className)){
        boolean diffOh = false;
        int oh = -1;
        Hashtable<String,Hashtable<String,LipidParameterSet>> resultsModIgnored = new Hashtable<String,Hashtable<String,LipidParameterSet>>();
        Hashtable<String,Hashtable<String,LipidParameterSet>> negatives = new Hashtable<String,Hashtable<String,LipidParameterSet>>();
        String anyValidRuleName = null;
        float maxDev = -1f;
        //TODO: here have to introduce separate RT processing according to the number of OHs!!!
        for (String mod : postProcessData.get(className).keySet()){
          try{
            if (RulesContainer.isRtPostprocessing(StaticUtils.getRuleName(className,mod))){
              anyValidRuleName = StaticUtils.getRuleName(className,mod);
              float aDev = -1f;
              if (RulesContainer.getRetentionTimeMaxDeviation(anyValidRuleName)!=null) aDev = new Float(RulesContainer.getRetentionTimeMaxDeviation(anyValidRuleName));
              if (aDev>0 && aDev>maxDev) maxDev = aDev;
            }
          }catch(Exception ex){}
          Hashtable<String,Hashtable<String,LipidParameterSet>> analytes = postProcessData.get(className).get(mod);
          Hashtable<String,Hashtable<String,LipidParameterSet>> negative = new Hashtable<String,Hashtable<String,LipidParameterSet>>();
          if (negativeExamples.containsKey(className) && negativeExamples.get(className).containsKey(mod)) negative = negativeExamples.get(className).get(mod);
          for (String analyteName : analytes.keySet()){
            Hashtable<String,LipidParameterSet> sameAnalyte = new Hashtable<String,LipidParameterSet>();
            if (resultsModIgnored.containsKey(analyteName)) sameAnalyte = resultsModIgnored.get(analyteName);
            for (String rt : analytes.get(analyteName).keySet()){
              if (sameAnalyte.containsKey(rt)){
                int count = 1;
                while (sameAnalyte.containsKey(rt+"_"+count)) count++;
                sameAnalyte.put(rt+"_"+String.valueOf(count), analytes.get(analyteName).get(rt));
              } else sameAnalyte.put(rt, analytes.get(analyteName).get(rt));
              
              //for oh
              if (oh==-1)
                oh = analytes.get(analyteName).get(rt).getOhNumber();
              else if (oh!=analytes.get(analyteName).get(rt).getOhNumber())
                diffOh = true;
            }
            resultsModIgnored.put(analyteName,sameAnalyte);
          }
          for (String analyteName : negative.keySet()){
            Hashtable<String,LipidParameterSet> sameAnalyte = new Hashtable<String,LipidParameterSet>();
            if (negatives.containsKey(analyteName)) sameAnalyte = negatives.get(analyteName);
            for (String rt : negative.get(analyteName).keySet()){
              if (sameAnalyte.containsKey(rt)){
                int count = 1;
                while (sameAnalyte.containsKey(rt+"_"+count)) count++;
                sameAnalyte.put(rt+"_"+String.valueOf(count), negative.get(analyteName).get(rt));
              } else sameAnalyte.put(rt, negative.get(analyteName).get(rt));
            }
            negatives.put(analyteName,sameAnalyte);
          }
        }
        float minDev = extractMinimumAcceptedDeviationValue(className,resultsModIgnored);
        try {
          resultsModIgnored = filterRetentionTimeSeries(className,anyValidRuleName,resultsModIgnored,negatives,maxDev,minDev,true,diffOh);
        }
        catch (LMException e) {
          System.out.println("Warning: "+className+" was not RT filtered: "+e.getMessage());
        }
        // store the results to the returning hash
        Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>> resultsClass = new Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>();
        for (String analyteName : resultsModIgnored.keySet()){
          Hashtable<String,Hashtable<String,LipidParameterSet>> resultsAnalyte = new Hashtable<String,Hashtable<String,LipidParameterSet>>();
          Hashtable<String,LipidParameterSet> sameAnalyte = resultsModIgnored.get(analyteName);
          for (String rt : sameAnalyte.keySet()){
            LipidParameterSet set = sameAnalyte.get(rt);
            String finalRt = new String(rt);
            if (finalRt.indexOf("_")!=-1) finalRt = finalRt.substring(0,finalRt.indexOf("_"));
            String mod = set.getModificationName();
            Hashtable<String,LipidParameterSet> resultsMod = new Hashtable<String,LipidParameterSet>();
            if (resultsAnalyte.containsKey(mod)) resultsMod = resultsAnalyte.get(mod);
            resultsMod.put(finalRt, set);
            resultsAnalyte.put(mod, resultsMod);
            resultsClass.put(analyteName, resultsAnalyte);
            results.put(className,resultsClass);
          }
        }
      }else{ 
        for (String mod : postProcessData.get(className).keySet()){
          Hashtable<String,Hashtable<String,LipidParameterSet>> negatives = new Hashtable<String,Hashtable<String,LipidParameterSet>>();
          if (negativeExamples.containsKey(className) && negativeExamples.get(className).containsKey(mod)) negatives = negativeExamples.get(className).get(mod);
          String ruleName = StaticUtils.getRuleName(className, mod);
          Hashtable<String,Hashtable<String,LipidParameterSet>> result = postProcessData.get(className).get(mod);
          //this is for oh
          boolean diffOh = false;
          int oh = -1;
          for (Hashtable<String,LipidParameterSet> sameAnalyte : result.values()){
            for (LipidParameterSet set : sameAnalyte.values()) {
              if (oh==-1)
                oh = set.getOhNumber();
              else if (oh!=set.getOhNumber())
                diffOh = true;
           }
          }
          try {
            float maxDev = -1f;
            if (RulesContainer.getRetentionTimeMaxDeviation(ruleName)!=null) maxDev = new Float(RulesContainer.getRetentionTimeMaxDeviation(ruleName));
            float minDev = extractMinimumAcceptedDeviationValue(className,postProcessData.get(className).get(mod));
            result = filterRetentionTimeSeries(className,ruleName,postProcessData.get(className).get(mod),negatives,maxDev,minDev,false,diffOh);
          }
          catch (LMException e) {
            System.out.println("Warning: "+className+"_"+mod+" was not RT filtered: "+e.getMessage());
          }

          // store the results to the returning hash
          Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>> resultsClass = new Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>();
          if (results.containsKey(className)) resultsClass = results.get(className);
          for (String analyteName : result.keySet()){
            Hashtable<String,Hashtable<String,LipidParameterSet>> resultsAnalyte = new Hashtable<String,Hashtable<String,LipidParameterSet>>();
            if (resultsClass.containsKey(analyteName)) resultsAnalyte = resultsClass.get(analyteName);
            resultsAnalyte.put(mod, result.get(analyteName));
            resultsClass.put(analyteName, resultsAnalyte);
            results.put(className,resultsClass);
          }
        }
      }
    }
    return results;
  }
  
  /**
   * orders data in a manner that is appropriate for LM predection
   * @param unprocessed the unprocessed data
   * @param ms2Removed data that is wrong according to MSn
   * @return the ordered data set
   */
  private Vector<Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>>> orderData (Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> unprocessed,
      Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> ms2Removed, Hashtable<String,Boolean> adductInsensitiveRtFilter){
    Vector<Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>>> dataOrdered = new Vector<Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>>>();
    Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> results = new Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>>();
    // first key class name, second key modification name, third key analyte name, last key retention time
    Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> postProcessData = new Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>>();
    Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> negativeExamples = new Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>>();
    
    // first, figure out which class-mod combination requires a post processing
    for (String className : unprocessed.keySet()){
      boolean adductInsensitive = false;
      if (adductInsensitiveRtFilter!=null && adductInsensitiveRtFilter.containsKey(className) && adductInsensitiveRtFilter.get(className)==true)
        adductInsensitive = true;
      boolean acceptAnyMod = false;
      for (String analyteName : unprocessed.get(className).keySet()){
        if (!adductInsensitive) continue;
        boolean modFound = false;
        for (String mod : unprocessed.get(className).get(analyteName).keySet()){
          try{
            if (RulesContainer.isRtPostprocessing(StaticUtils.getRuleName(className,mod))){
              modFound = true;
              break;
            }
          }catch(Exception ex){}
        }
        if (modFound){
          acceptAnyMod = true;
          break;
        }
      }
      Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>> postProcessOfClass = new Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>();
      for (String analyteName : unprocessed.get(className).keySet()){
        for (String mod : unprocessed.get(className).get(analyteName).keySet()){
          try{
            if (acceptAnyMod || RulesContainer.isRtPostprocessing(StaticUtils.getRuleName(className,mod))){
              Hashtable<String,Hashtable<String,LipidParameterSet>> postProcessOfMod = new Hashtable<String,Hashtable<String,LipidParameterSet>>();
              if (postProcessOfClass.containsKey(mod)) postProcessOfMod = postProcessOfClass.get(mod);
              postProcessOfMod.put(analyteName, unprocessed.get(className).get(analyteName).get(mod));
              postProcessOfClass.put(mod, postProcessOfMod);
              
              
            }
          }catch(Exception ex){}
          // if there is no post processing for this class and modification - the results can be put to the results hash
          if (!postProcessOfClass.containsKey(mod)){
            Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>> resultsClass = new Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>();
            if (results.containsKey(className)) resultsClass = results.get(className);
            Hashtable<String,Hashtable<String,LipidParameterSet>> resultsAnalyte = new Hashtable<String,Hashtable<String,LipidParameterSet>>();
            if (resultsClass.containsKey(analyteName)) resultsAnalyte = resultsClass.get(analyteName);
            resultsAnalyte.put(mod, unprocessed.get(className).get(analyteName).get(mod));
            resultsClass.put(analyteName, resultsAnalyte);
            results.put(className,resultsClass);
          }
        }
      }
      if (postProcessOfClass.size()>0){
        postProcessData.put(className, postProcessOfClass);
        // if we have data to post process for a class and mod - we have to extract the negative examples for them
        Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>> negativeClass = new Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>();
        if (ms2Removed.containsKey(className)){
          for (String analyteName : ms2Removed.get(className).keySet()){
            for (String mod: ms2Removed.get(className).get(analyteName).keySet()){
              if (!postProcessOfClass.containsKey(mod)) continue;
              Hashtable<String,Hashtable<String,LipidParameterSet>> negativeMod = new Hashtable<String,Hashtable<String,LipidParameterSet>>();
              if (negativeClass.containsKey(mod)) negativeMod = negativeClass.get(mod);
              negativeMod.put(analyteName, ms2Removed.get(className).get(analyteName).get(mod));
              negativeClass.put(mod, negativeMod);
            }
          }
        }
        negativeExamples.put(className, negativeClass);
      } else if (!results.containsKey(className)){
        results.put(className, unprocessed.get(className));
      }
    }
    dataOrdered.add(results);
    dataOrdered.add(postProcessData);
    dataOrdered.add(negativeExamples);
    return dataOrdered;
  }
  
  /**
   * this method extracts the MSn hits to use first, then, it starts the iterative
   * Levenberg-Marquardt filtering process
   * @param className required for adduct insensitive filtering (only the ones where RetentionTimePostprocessing=true should be used for the model)
   * @param ruleName required to extract the number of C atoms and double bonds
   * @param unprocessed the unprocessed analytes 
   * @param negatives hits removed by MSn - useable for counter curve
   * @param maxDev maximally allowed deviation from the model
   * @param adductInsensitive was adduct insensitive processing chosen
   * @param respectOh has the OH to be fitted too
   * @return hits that pass the filter
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the library is not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * @throws LMException exception if model adaption is not possible
   */
  private Hashtable<String,Hashtable<String,LipidParameterSet>> filterRetentionTimeSeries(String className, String ruleName, Hashtable<String,Hashtable<String,LipidParameterSet>> unprocessed,
      Hashtable<String,Hashtable<String,LipidParameterSet>> negatives, float maxDev, float minDev, boolean adductInsensitive, boolean respectOh) throws RulesException, NoRuleException, IOException, SpectrummillParserException, LMException{
    // find the MSn identifications and group them according to their # C atoms and # double bonds
    Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> paramsOrdered = new Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>();
    Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> unprocessedOrdered = new Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>();
    @SuppressWarnings("rawtypes")
    Vector ranges = groupAccordingToCAtomsAndDoubleBonds(unprocessed,paramsOrdered,unprocessedOrdered,className,ruleName);
    RangeInteger cAtomsRange = (RangeInteger)ranges.get(0);
    @SuppressWarnings("unchecked")
    Hashtable<Integer,RangeInteger> dbsRanges = (Hashtable<Integer,RangeInteger>)ranges.get(1);
    Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> negativesOrdered = new Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>();
    @SuppressWarnings({ "rawtypes", "unused" })
    Vector rangesNeg = groupAccordingToCAtomsAndDoubleBonds(negatives,new Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>(),negativesOrdered,className,ruleName);
    float tolerance = 4f;
    
    @SuppressWarnings("rawtypes")
    Vector results = doIterativeLMOptimization(unprocessedOrdered,paramsOrdered,cAtomsRange,dbsRanges,negativesOrdered,tolerance,maxDev,respectOh);
    LevenbergMarquardtOptimizer optimizer = (LevenbergMarquardtOptimizer)results.get(0);
    if (adductInsensitive) {
      predictedModels_.put(className, optimizer);
      respectOhs_.put(className, respectOh);
    }else {
      predictedModels_.put(ruleName, optimizer);
      respectOhs_.put(ruleName, respectOh);
    }
    LevenbergMarquardtOptimizer counterModel = (LevenbergMarquardtOptimizer)results.get(1);
    @SuppressWarnings("unchecked")
    Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> foundNegatives = (Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>)results.get(2);

    Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> newForModel = filterValidHitsBasedOnModel(optimizer,
        unprocessedOrdered,cAtomsRange,dbsRanges,tolerance,maxDev,minDev,counterModel,foundNegatives,respectOh);
    
    Hashtable<String,Hashtable<String,LipidParameterSet>> filteredValues = new Hashtable<String,Hashtable<String,LipidParameterSet>>();
    Pattern cAtomsPattern =  Pattern.compile(RulesContainer.getCAtomsFromNamePattern(ruleName));
    Pattern dbsPattern =  Pattern.compile(RulesContainer.getDoubleBondsFromNamePattern(ruleName));
    for (String analyteName : unprocessed.keySet()){
      Hashtable<String,LipidParameterSet> sameAnalyte = unprocessed.get(analyteName);
      Matcher cAtomsMatcher = cAtomsPattern.matcher(analyteName.split(LipidomicsConstants.CHAIN_MOD_SEPARATOR)[0]);
      if (!cAtomsMatcher.matches()) throw new RulesException("The analyte "+analyteName+" does not match the "+FragRuleParser.GENERAL_CATOMS_PARSE+" pattern \""+RulesContainer.getCAtomsFromNamePattern(ruleName)+"\" of the class "+ruleName+"!");
      int cAtoms = Integer.parseInt(cAtomsMatcher.group(1));
      Matcher dbsMatcher = dbsPattern.matcher(analyteName.split(LipidomicsConstants.CHAIN_MOD_SEPARATOR)[0]);
      if (!dbsMatcher.matches()) throw new RulesException("The analyte "+analyteName+" does not match the "+FragRuleParser.GENERAL_DBOND_PARSE+" pattern \""+RulesContainer.getDoubleBondsFromNamePattern(ruleName)+"\" of the class "+ruleName+"!");
      int dbs = Integer.parseInt(dbsMatcher.group(1));
      if (!newForModel.containsKey(cAtoms) || !newForModel.get(cAtoms).containsKey(dbs)) continue;
      Hashtable<String,LipidParameterSet> filteredRts = newForModel.get(cAtoms).get(dbs);
      Hashtable<String,LipidParameterSet> acceptedRts = new Hashtable<String,LipidParameterSet>();
      for (String rt : sameAnalyte.keySet()){
        if (filteredRts.containsKey(rt)) acceptedRts.put(rt, sameAnalyte.get(rt));
      }
      if (acceptedRts.size()>0) filteredValues.put(analyteName, acceptedRts);
    }
//    System.out.println("Params: ");
//    FloatMatrix params = optimizer.getResultParams();
//    for (int i=0;i!=params.m;i++){
//      System.out.println(params.A[i][0]);
//    }
//    for (Integer cAtoms : newForModel.keySet()){
//      for (Integer dbs : newForModel.get(cAtoms).keySet()){
//        Hashtable<String,LipidParameterSet> sameRt = newForModel.get(cAtoms).get(dbs);
//        for (String rt : sameRt.keySet()){
//          System.out.println("MS2-Hit: "+cAtoms+":"+dbs+"_"+rt);
//        }
//      }      
//    }  

    return filteredValues;
  }
  
  /**
   * performs iterative Levenberg-Marquardt curve optimization - iterative because the filter takes
   * iteratively neighboring carbon atoms and double bonds into account
   * @param unprocessed data to be filtered
   * @param forModel data used for model optimization
   * @param cAtomsMaxRange the maximum range of carbon atoms for which the model shall be optimized
   * @param dbsMaxRanges the maximum range of double bonds for which the model shall be optimized
   * @param negatives MS/MS hits that do not belong to this class can be used to generate a counter model to remove hits
   * @param tolerance the multiplication factor for the mean deviation
   * @param maxDev the maximum RT deviation allowed
   * @param respectOh has the OH to be fitted too
   * @return the optimized model in form of Levenberg-Marquardt object
   * @throws LMException exception if model adaption is not possible
   */
  @SuppressWarnings("rawtypes")
  private Vector doIterativeLMOptimization(Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> unprocessed,
      Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> forModel, RangeInteger cAtomsMaxRange, Hashtable<Integer,RangeInteger> dbsMaxRanges,
      Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> negatives, float tolerance, float maxDev,
      boolean respectOh) throws LMException{
    return doIterativeLMOptimization(unprocessed, forModel, cAtomsMaxRange, dbsMaxRanges, null, null, null,negatives,
        new Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>(),tolerance,maxDev,0,respectOh);
  }
  
  @SuppressWarnings({ "unchecked", "rawtypes" })
  private Vector doIterativeLMOptimization(Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> unprocessed,
      Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> forModel, RangeInteger cAtomsMaxRange, Hashtable<Integer,RangeInteger> dbsMaxRanges,
      RangeInteger cAtomsRange, Hashtable<Integer,RangeInteger> dbsRanges, LevenbergMarquardtOptimizer prevOptimizer,
      Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> negatives, 
      Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> remainingNegatives, float tolerance, float maxDev, int count,
      boolean respectOh) throws LMException{
    // the first step is to optimize the model
//    FloatMatrix prevParams = null;
//    float prevMeanDeviation = Float.MAX_VALUE;
    if (prevOptimizer!=null){
//      prevParams = prevOptimizer.getResultParams();
//      prevMeanDeviation = prevOptimizer.getMeanDeviation();
    }
    LevenbergMarquardtOptimizer optimizer = optimizeLMModel(forModel, null,respectOh);
//    if (optimizer.getMeanDeviation()>prevMeanDeviation) optimizer.setMaxDeviation(prevMeanDeviation);
    //System.out.println(optimizer.getResultChiSqr()+";"+optimizer.getMeanDeviation()+";"+optimizer.getResultLambda());
    
    // second, figure out which C and double bond range should be filtered by the calculated model
    Vector currentRanges = getCurrentRanges(forModel, cAtomsRange, dbsRanges);
    RangeInteger currentCAtomsRange = (RangeInteger)currentRanges.get(0);
    Hashtable<Integer,RangeInteger> currentDbsRanges = (Hashtable<Integer,RangeInteger>)currentRanges.get(1);
    
    
    LevenbergMarquardtOptimizer counterModel = null;
    Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> returnedNegatives = new Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>();
    Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> unusedNegatives = new Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>();
    try{
      Vector result = checkForLMCounterModel(optimizer,negatives,remainingNegatives,currentCAtomsRange,currentDbsRanges,tolerance,respectOh);
       if (result!=null){
         counterModel = (LevenbergMarquardtOptimizer)result.get(0);
         returnedNegatives = (Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>)result.get(1);
         unusedNegatives = (Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>)result.get(2);
       }
    } catch (LMException lmx){
      System.out.println("Warning: calculation of counter model was not possible! Reason: "+lmx.getMessage());
    }
    if (counterModel!=null && count==0){
      Vector<LevenbergMarquardtOptimizer> optimizers = evaluateFilterAgainBasedOnCounterModel(forModel,returnedNegatives,optimizer,counterModel,respectOh);
      

      optimizer = optimizers.get(0);
      counterModel = optimizers.get(1);

    }  
    
    
    Vector ranges = proposeLMFilterRanges(unprocessed, cAtomsMaxRange, dbsMaxRanges, currentCAtomsRange, currentDbsRanges);
    RangeInteger newCAtomsRange = (RangeInteger)ranges.get(0);
    Hashtable<Integer,RangeInteger> newDbsRanges = (Hashtable<Integer,RangeInteger>) ranges.get(1);
    
    if (counterModel!=null){
      Vector result = checkForLMCounterModel(optimizer,returnedNegatives,unusedNegatives,newCAtomsRange,newDbsRanges,tolerance,respectOh);
      if (result!=null && result.size()>0 && result.get(0)!=null){
        counterModel = (LevenbergMarquardtOptimizer)result.get(0);
        returnedNegatives = (Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>)result.get(1);
        unusedNegatives = (Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>)result.get(2);
      }
    }
        
    // optimizer covers all ranges -> optimization finished   
////    if (allRangesCovered(currentDbsRanges,dbsMaxRanges)){
      
//      System.out.println("Params: ");
//      FloatMatrix params = optimizer.getResultParams();
//      for (int i=0;i!=params.m;i++){
//        System.out.println(params.A[i][0]);
//      }
      Vector results = new Vector();
      results.add(optimizer);
      results.add(counterModel);
      results.add(returnedNegatives);
      
      return results;
////    } else{
      // apply the filter and fetch hits for a new model - filter moves to next C atoms and double bond range
////      Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> newForModel = filterValidHitsBasedOnModel(optimizer,
////          unprocessed,newCAtomsRange,newDbsRanges,tolerance,maxDev,counterModel,returnedNegatives);
////      return doIterativeLMOptimization(unprocessed, newForModel, cAtomsMaxRange, dbsMaxRanges, newCAtomsRange,
////          newDbsRanges, optimizer,returnedNegatives,unusedNegatives,tolerance,maxDev,count+1);
////    }    
  }
  
  /**
   * optimizes the Levenberg-Marquardt model - data is transformed in adequate structure, optimizer initialized and model fitted
   * @param forModel data to be used for the Levenberg-Marquardt model 
   * @param prevParams params of a previous Levenberg-Marquardt model (to speed up optimization)
   * @param respectOh has the OH to be fitted too
   * @return Levenberg-Marquardt optimizer class including the fitted model
   * @throws LMException exception if model adaption is not possible
   */
  private LevenbergMarquardtOptimizer optimizeLMModel(Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> forModel, FloatMatrix prevParams, boolean respectOh) throws LMException{
    int dataSize = 0;
    for (Hashtable<Integer,Hashtable<String,LipidParameterSet>> values1 : forModel.values()){
      for (Hashtable<String,LipidParameterSet> values2 : values1.values()){
        for (@SuppressWarnings("unused") LipidParameterSet set : values2.values()) dataSize++;
      }
    }
    float[][] values;
    if (respectOh)
      values = new float[dataSize][3];
    else
      values = new float[dataSize][2];
    float[] rtValues = new float[dataSize];
    int count = 0;
    for (Integer cAtoms : forModel.keySet()){
      Hashtable<Integer,Hashtable<String,LipidParameterSet>> sameCAtoms = forModel.get(cAtoms);
      for (Integer dbs : sameCAtoms.keySet()){
        Hashtable<String,LipidParameterSet> sameCAndDbs = sameCAtoms.get(dbs);
        for (String rt : sameCAndDbs.keySet()){
          values[count][0] = cAtoms;
          values[count][1] = dbs;
          if (respectOh)
            values[count][2] = sameCAndDbs.get(rt).getOhNumber();
          rtValues[count] = getRtFromRtString(rt);
          count++; 
        }
      }
    }
    LevenbergMarquardtOptimizer lmOptimizer = null;
    if (respectOh) {
      try{
        lmOptimizer = new LMAsymptVarDecayThreeVariables(values,rtValues,prevParams);
        lmOptimizer.fit();
        //System.out.println("I used the three variables");
      } catch (Exception ex){
        try{
          //System.out.println("Warning: trying fit with 1/x");
          lmOptimizer = new LMAsymptDecayThreeVariables(values,rtValues,prevParams);
          lmOptimizer.fit();
        } catch (Exception ex2){
          //System.out.println("Warning: trying fit with log(x)");
          try{
            lmOptimizer = new LMLogDecayThreeVariables(values,rtValues,prevParams);
            lmOptimizer.fit();
          } catch (Exception ex3){
            throw new LMException("The model could not be fitted!");
          }
        }
      }      
    } else {
      try{
        lmOptimizer = new LMAsymptVarDecayTwoVariables(values,rtValues,prevParams);
        //LevenbergMarquardtOptimizer lmOptimizer = new LMLogDecayTwoVariables(values,rtValues,prevParams);
        //LevenbergMarquardtOptimizer lmOptimizer = new LMEulerTwoVariables(values,rtValues,prevParams);
        //LevenbergMarquardtOptimizer lmOptimizer = new LMQuadraticTwoVariables(values,rtValues,prevParams);
        lmOptimizer.fit();
      } catch (Exception ex){
        try{
          //System.out.println("Warning: trying fit with 1/x");
          lmOptimizer = new LMAsymptDecayTwoVariables(values,rtValues,prevParams);
          lmOptimizer.fit();
        } catch (Exception ex2){
          //System.out.println("Warning: trying fit with log(x)");
          try{
            lmOptimizer = new LMLogDecayTwoVariables(values,rtValues,prevParams);
            lmOptimizer.fit();
          } catch (Exception ex3){
            throw new LMException("The model could not be fitted!");
          }
        }
      }
    }
    return lmOptimizer;
  }
  
  /**
   * extracts currently used ranges for carbon atoms and double bonds for the model
   * @param forModel data that was used for the model
   * @param cAtomsRange predefined carbon atoms range (set to null if you want to extract from the data)
   * @param dbsRanges predefined double bond range (set to null if you want to extract from the data)
   * @return currently used ranges for carbon atoms and double bonds for the model
   */
  @SuppressWarnings({ "rawtypes", "unchecked" })
  private Vector getCurrentRanges(Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> forModel,
      RangeInteger cAtomsRange, Hashtable<Integer,RangeInteger> dbsRanges){
    Vector results = new Vector();
    RangeInteger currC = null;
    Hashtable<Integer,RangeInteger> currDbs = null;
    if (cAtomsRange == null || dbsRanges ==null){
      currC = new RangeInteger(Integer.MAX_VALUE,0);
      currDbs = new Hashtable<Integer,RangeInteger>();
      for (Integer cAtoms : forModel.keySet()){
        if (cAtoms<currC.getStart()) currC = new RangeInteger(cAtoms,currC.getStop());
        if (cAtoms>currC.getStop()) currC = new RangeInteger(currC.getStart(),cAtoms);
        for (Integer dbs : forModel.get(cAtoms).keySet()){
          RangeInteger dbsRange = new RangeInteger(Integer.MAX_VALUE,0);
          if (currDbs.containsKey(cAtoms)) dbsRange = currDbs.get(cAtoms);
          if (dbs<dbsRange.getStart()) dbsRange = new RangeInteger(dbs,dbsRange.getStop());
          if (dbs>dbsRange.getStop()) dbsRange = new RangeInteger(dbsRange.getStart(),dbs);
          currDbs.put(cAtoms, dbsRange);
        }
      }
    } else {
      currC = cAtomsRange;
      currDbs = dbsRanges;
    }
    results.add(currC);
    results.add(currDbs);
    return results;
  }
  
  /**
   * proposes carbon atoms and double bond ranges to which the Levenberg-Marquardt filter shall be applied
   * @param unprocessed available data for filtering
   * @param cAtomsMaxRange the maximum range of carbon atoms for which the model shall be optimized
   * @param dbsMaxRanges the maximum range of double bonds for which the model shall be optimized
   * @param cAtomsRange the range of carbon atoms that is used for Levenberg-Marquardt filter optimization
   * @param dbsRanges the range of double bonds that is used for Levenberg-Marquardt filter optimization
   * @return a Vector containing the proposed carbon atoms range at first position, and the proposed double bonds range at second position
   */
  @SuppressWarnings({ "rawtypes", "unchecked" })
  private Vector proposeLMFilterRanges(Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> unprocessed,
      RangeInteger cAtomsMaxRange, Hashtable<Integer,RangeInteger> dbsMaxRanges,
      RangeInteger cAtomsRange, Hashtable<Integer,RangeInteger> dbsRanges){
    
    // find the closest C atom distance and dbs of a unprocessed result to the current ranges
    boolean oneMissing = false;
    boolean allCRangesCovered = true;
    int minCDistance = Integer.MAX_VALUE;
    int minDbsDistance = Integer.MAX_VALUE;
    for (Integer cAt : unprocessed.keySet()){   
      if (dbsRanges.containsKey(cAt)){
        RangeInteger dbRange = dbsRanges.get(cAt);
        Hashtable<Integer,Hashtable<String,LipidParameterSet>> sameC = unprocessed.get(cAt);
        for (Integer dbs : sameC.keySet()){
          if (dbRange.insideRange(dbs)) continue;
          oneMissing = true;
          int dist = dbs-dbRange.getStop();
          if (dbs<dbRange.getStart()) dist = dbRange.getStart()-dbs;
          if (dist<minDbsDistance) minDbsDistance = dist;
        }
      // if the no analyte with this C atoms count is in the list -> determine the min C distance        
      } else {
        allCRangesCovered = false;
        oneMissing = true;
        for (Integer cFound : dbsRanges.keySet()){
          int dist = Math.abs(cAt-cFound);
          if (dist<minCDistance) minCDistance = dist;
        }
      }
    }
    // now propose the new ranges
    Vector results = new Vector();
    if (oneMissing){
      Hashtable<Integer,RangeInteger> proposedDbs = new Hashtable<Integer,RangeInteger>();
      List<Integer> cAtomsList = new ArrayList<Integer>(unprocessed.keySet());
      Collections.sort(cAtomsList);
      for (int i=0; i!=cAtomsList.size(); i++){
        Integer cAtoms = cAtomsList.get(i);
        Hashtable<Integer,Hashtable<String,LipidParameterSet>> sameC = unprocessed.get(cAtoms);
        RangeInteger dbsRange = null;
        if (dbsRanges.containsKey(cAtoms)){
          dbsRange = dbsRanges.get(cAtoms);
          if (sameC.containsKey(dbsRange.getStart()-minDbsDistance)){
            dbsRange = new RangeInteger(dbsRange.getStart()-minDbsDistance,dbsRange.getStop());
          }
          if (sameC.containsKey(dbsRange.getStop()+minDbsDistance)){
            dbsRange = new RangeInteger(dbsRange.getStart(),dbsRange.getStop()+minDbsDistance);
          }
        }
        RangeInteger maxDbsRange = dbsMaxRanges.get(cAtoms);
        // is there a found hit with less C atoms
        if (i!=0 && (dbsRanges.containsKey(cAtoms-minCDistance) || allCRangesCovered)){
          RangeInteger lessCs = null;
          if (allCRangesCovered) lessCs = dbsRanges.get(cAtomsList.get(i-1));
          else lessCs = dbsRanges.get(cAtoms-minCDistance);
          if (dbsRange == null) dbsRange = lessCs;
          else dbsRange.extendToOtherRanges(dbsRange,maxDbsRange);
        }
        // is there a found hit with more C atoms
        if (i!=(cAtomsList.size()-1) && (dbsRanges.containsKey(cAtoms+minCDistance) || allCRangesCovered)){
          RangeInteger moreCs = null;
          if (allCRangesCovered) moreCs =  dbsRanges.get(cAtomsList.get(i+1));
          else moreCs =  dbsRanges.get(cAtoms+minCDistance);
          if (dbsRange == null) dbsRange = moreCs;
          else dbsRange.extendToOtherRanges(dbsRange,maxDbsRange);
        }        
        if (dbsRange!=null) proposedDbs.put(cAtoms, dbsRange);
      }
      int lowestCs = Integer.MAX_VALUE;
      int highestCs = 0;
      for (int cAtoms : proposedDbs.keySet()){
        if (cAtoms<lowestCs) lowestCs = cAtoms;
        if (cAtoms>highestCs) highestCs = cAtoms;
      }
      RangeInteger proposedCs = new RangeInteger(lowestCs,highestCs);
      results.add(proposedCs);
      results.add(proposedDbs);
    } else {

      results.add(cAtomsMaxRange);
      results.add(dbsMaxRanges);
    }    
    return results;
  }
  
  /**
   * filters the analytes based on the optimized model within the proposed ranges
   * @param model mathematical model returned from the Levenberg-Marquardt optimization process
   * @param unprocessed unfiltered data
   * @param cAtomsRange the carbon atoms range to which the filter shall be applied
   * @param dbsRanges the carbon atoms range to which the filter shall be applied
   * @param tolerance the multiplication factor for the mean deviation
   * @param maxDev the maximum RT deviation allowed
   * @param counterModel the model fitted based on the negative hits
   * @param negatives MS/MS hits that do not belong to this class can be used to generate a counter model to remove hits
   * @return data passing the filter
   * @throws LMException exception if model adaption is not possible
   */
  private Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> filterValidHitsBasedOnModel(LevenbergMarquardtOptimizer model,
      Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> unprocessed,
      RangeInteger cAtomsRange, Hashtable<Integer,RangeInteger> dbsRanges, float tolerance, float maxDev, float minDev,
      LevenbergMarquardtOptimizer counterModel, Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> negatives,
      boolean respectOh) throws LMException {
    Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> forModel = new Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>();
//    System.out.println("C-atoms range: "+cAtomsRange.getStart()+"-"+cAtomsRange.getStop());
//    System.out.println("MeanDev: "+model.getMeanDeviation());
    for (Integer cAtoms : unprocessed.keySet()){
      if (!dbsRanges.containsKey(cAtoms)) continue;
      Hashtable<Integer,Hashtable<String,LipidParameterSet>> sameC = unprocessed.get(cAtoms);
      Hashtable<Integer,Hashtable<String,LipidParameterSet>> sameCModel = new Hashtable<Integer,Hashtable<String,LipidParameterSet>>();
      Hashtable<Integer,Hashtable<String,LipidParameterSet>> negSameC = new Hashtable<Integer,Hashtable<String,LipidParameterSet>>();
      if (negatives.containsKey(cAtoms)) negSameC = negatives.get(cAtoms);
      RangeInteger dbsRange = dbsRanges.get(cAtoms);
//      System.out.println(cAtoms+" Dbs-Range: "+dbsRange.getStart()+" ; "+dbsRange.getStop());
      for (Integer dbs : sameC.keySet()){
        if (!dbsRange.insideRange(dbs)) continue;
        Hashtable<String,LipidParameterSet> sameDbs = sameC.get(dbs);
        Hashtable<String,LipidParameterSet> sameDbsModel = new Hashtable<String,LipidParameterSet>();
        Hashtable<String,LipidParameterSet> negSameDbs = new Hashtable<String,LipidParameterSet>();
        if (negSameC.containsKey(dbs)) negSameDbs = negSameC.get(dbs);
        for (String rtString : sameDbs.keySet()){
          if (negSameDbs.containsKey(rtString)){
            continue;
          }
          float rt = getRtFromRtString(rtString);
          float[] cAndDbs = new float[2];
          if (respectOh)
            cAndDbs = new float[3];
          cAndDbs[0] = cAtoms;
          cAndDbs[1] = dbs;
          if (respectOh)
            cAndDbs[2] = sameDbs.get(rtString).getOhNumber();
          float fitRt = model.calculateFitValue(cAndDbs);
          float counterRt = Float.NaN;
          if (counterModel!=null) counterRt = counterModel.calculateFitValue(cAndDbs);
          float allowedRange = model.getMeanDeviation()*tolerance;
          if (maxDev>0 && allowedRange>maxDev) allowedRange = maxDev;
          if (minDev>allowedRange) allowedRange = minDev;
          if ((fitRt-allowedRange)<rt && rt<(fitRt+allowedRange)){
            // checks if the hit is closer to the counter model
            if (fallsOnCounterModelSide(rt,fitRt,counterRt)){
              negSameDbs.put(rtString, sameDbs.get(rtString));
              //System.out.println("Is a FP: "+cAtoms+":"+dbs+" "+rt+" "+fitRt+";"+counterRt);
            }else
              sameDbsModel.put(rtString, sameDbs.get(rtString));
          }
        }
        if (sameDbsModel.size()>0) sameCModel.put(dbs, sameDbsModel);
        if (negSameDbs.size()>0) negSameC.put(dbs, negSameDbs);
      }
      if (sameCModel.size()>0) forModel.put(cAtoms, sameCModel);
      if (negSameC.size()>0) negatives.put(cAtoms, negSameC);
    }
    return forModel;
  }
  
  /**
   * checks if the algorithm reached the max range - if so, stop at this iteration
   * @param currentDbsRanges the current double bonds ranges to which the filter shall be applied
   * @param dbsMaxRanges the maximum double bonds ranges to which the filter shall be applied
   * @return true if both ranges are the same
   */
//  private boolean allRangesCovered(Hashtable<Integer,RangeInteger> currentDbsRanges, Hashtable<Integer,RangeInteger> dbsMaxRanges){
//    boolean coversAllRanges = true;
//    for (Integer cAtoms : dbsMaxRanges.keySet()){
//      if (!currentDbsRanges.containsKey(cAtoms)){
//        coversAllRanges = false;
//        break;
//      }
//      RangeInteger currentDbs = currentDbsRanges.get(cAtoms);
//      RangeInteger dbsMaxRange = dbsMaxRanges.get(cAtoms);
//      if (currentDbs.getStart()>dbsMaxRange.getStart() || currentDbs.getStop()<dbsMaxRange.getStop()){
//        coversAllRanges = false;
//        break;        
//      }
//    }
//    return coversAllRanges;
//  }
  
  @SuppressWarnings({ "rawtypes", "unchecked" })
  private Vector groupAccordingToCAtomsAndDoubleBonds(Hashtable<String,Hashtable<String,LipidParameterSet>> input, 
      Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> msnOrdered,
      Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> allOrdered,
      String className, String ruleName) throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    Pattern cAtomsPattern =  Pattern.compile(RulesContainer.getCAtomsFromNamePattern(ruleName));
    Pattern dbsPattern =  Pattern.compile(RulesContainer.getDoubleBondsFromNamePattern(ruleName));
    RangeInteger cAtomsRange = new RangeInteger(Integer.MAX_VALUE,0);
    Hashtable<Integer,RangeInteger> dbsRanges = new Hashtable<Integer,RangeInteger>();

    for (String analyteName : input.keySet()){
    	if (analyteName.contains(LipidomicsConstants.CHAIN_MOD_SEPARATOR))
    		System.out.println("ox");
      Matcher cAtomsMatcher = cAtomsPattern.matcher(analyteName.split(LipidomicsConstants.CHAIN_MOD_SEPARATOR)[0]);
      if (!cAtomsMatcher.matches()) throw new RulesException("The analyte "+analyteName+" does not match the "+FragRuleParser.GENERAL_CATOMS_PARSE+" pattern \""+RulesContainer.getCAtomsFromNamePattern(ruleName)+"\" of the class "+ruleName+"!");
      int cAtoms = Integer.parseInt(cAtomsMatcher.group(1));
      Matcher dbsMatcher = dbsPattern.matcher(analyteName.split(LipidomicsConstants.CHAIN_MOD_SEPARATOR)[0]);
      if (!dbsMatcher.matches()) throw new RulesException("The analyte "+analyteName+" does not match the "+FragRuleParser.GENERAL_DBOND_PARSE+" pattern \""+RulesContainer.getDoubleBondsFromNamePattern(ruleName)+"\" of the class "+ruleName+"!");
      int dbs = Integer.parseInt(dbsMatcher.group(1));
      if (cAtoms<cAtomsRange.getStart()) cAtomsRange = new RangeInteger(cAtoms,cAtomsRange.getStop());
      if (cAtoms>cAtomsRange.getStop()) cAtomsRange = new RangeInteger(cAtomsRange.getStart(),cAtoms);
      RangeInteger dbsRange = new RangeInteger(Integer.MAX_VALUE,0);
      if (dbsRanges.containsKey(cAtoms)) dbsRange = dbsRanges.get(cAtoms);
      if (dbs<dbsRange.getStart()) dbsRange = new RangeInteger(dbs,dbsRange.getStop());
      if (dbs>dbsRange.getStop()) dbsRange = new RangeInteger(dbsRange.getStart(),dbs);
      dbsRanges.put(cAtoms, dbsRange);
      
      Hashtable<String,LipidParameterSet> hitsOfAnalyte = input.get(analyteName);
      for (String rt : hitsOfAnalyte.keySet()){
        LipidParameterSet params = hitsOfAnalyte.get(rt);
        if (params.isSuitableForRtProcessingHit(className)){
          Hashtable<Integer,Hashtable<String,LipidParameterSet>> paramsSameC = new Hashtable<Integer,Hashtable<String,LipidParameterSet>>();
          if (msnOrdered.containsKey(cAtoms)) paramsSameC = msnOrdered.get(cAtoms);
          Hashtable<String,LipidParameterSet> paramsSameDbs = new  Hashtable<String,LipidParameterSet>();
          if (paramsSameC.containsKey(dbs)) paramsSameDbs = paramsSameC.get(dbs);
          paramsSameDbs.put(rt, params);
          paramsSameC.put(dbs,paramsSameDbs);
          msnOrdered.put(cAtoms, paramsSameC);
        }
        Hashtable<Integer,Hashtable<String,LipidParameterSet>> unprocSameC = new Hashtable<Integer,Hashtable<String,LipidParameterSet>>();
        if (allOrdered.containsKey(cAtoms)) unprocSameC = allOrdered.get(cAtoms);
        Hashtable<String,LipidParameterSet> unprocSameDbs = new  Hashtable<String,LipidParameterSet>();
        if (unprocSameC.containsKey(dbs)) unprocSameDbs = unprocSameC.get(dbs);
        unprocSameDbs.put(rt, params);
        unprocSameC.put(dbs,unprocSameDbs);
        allOrdered.put(cAtoms, unprocSameC);

      }
    }
    Vector results = new Vector();
    results.add(cAtomsRange);
    results.add(dbsRanges);
    return results;
  }
  
  @SuppressWarnings({ "rawtypes", "unchecked" })
  private Vector checkForLMCounterModel(LevenbergMarquardtOptimizer optimizer, Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> negatives,
      Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> remainingNegatives,
      RangeInteger cAtomsRange, Hashtable<Integer,RangeInteger> dbsRanges, float tolerance, boolean respectOh) throws LMException{
    Vector result = null;
    Vector<Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>> results = filterAdequateNegatives(optimizer,negatives,remainingNegatives,cAtomsRange,dbsRanges,tolerance);
    if (results!=null){
      try{
        Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> usables = results.get(0);
        Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> notUsed = results.get(1);
        LevenbergMarquardtOptimizer counterModel = optimizeLMModel(usables, null, respectOh);
        result = new Vector();
        result.add(counterModel);
        result.add(usables);
        result.add(notUsed);
      }catch(Exception ex){ }
    }
    return result;
  }
  
  private Vector<Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>> filterAdequateNegatives(LevenbergMarquardtOptimizer optimizer, Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> negatives,
      Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> remainingNegatives,
      RangeInteger cAtomsRange, Hashtable<Integer,RangeInteger> dbsRanges, float tolerance) throws LMException{
    Vector<Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>> results = new Vector<Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>>();
    Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> negs = new Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>();
    Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> unusedNegs = new Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>();
    Hashtable<Integer,Hashtable<Integer,Hashtable<String,Float>>> diffs = new Hashtable<Integer,Hashtable<Integer,Hashtable<String,Float>>>();
    
    // merge negatives and remaining negatives
    for (Integer cAtoms : remainingNegatives.keySet()){
      Hashtable<Integer,Hashtable<String,LipidParameterSet>> sameC = new Hashtable<Integer,Hashtable<String,LipidParameterSet>>();
      if (negatives.containsKey(cAtoms)) sameC = negatives.get(cAtoms);
      for (Integer dbs : remainingNegatives.get(cAtoms).keySet()){
        Hashtable<String,LipidParameterSet> sameDbs = new Hashtable<String,LipidParameterSet>();
        if (sameC.containsKey(dbs)) sameDbs = sameC.get(dbs);
        for (String rt : remainingNegatives.get(cAtoms).get(dbs).keySet()){
          sameDbs.put(rt, remainingNegatives.get(cAtoms).get(dbs).get(rt));
        }
        
        sameC.put(dbs, sameDbs);
      }
      negatives.put(cAtoms, sameC);
    }
    
    // check for negatives inside the tolerance and if it has a positive or negative deviation
    int countPos = 0;
    int countNeg = 0;
    for (Integer cAt : negatives.keySet()){
      Hashtable<Integer,Hashtable<String,LipidParameterSet>> sameC = negatives.get(cAt);
      if (!cAtomsRange.insideRange(cAt) || !dbsRanges.containsKey(cAt)){
        unusedNegs.put(cAt, sameC);
        continue;
      }
      Hashtable<Integer,Hashtable<String,Float>> sameCDiffs = new Hashtable<Integer,Hashtable<String,Float>>();
      Hashtable<Integer,Hashtable<String,LipidParameterSet>> unusedSameC = new Hashtable<Integer,Hashtable<String,LipidParameterSet>>();
      RangeInteger dbsRange = dbsRanges.get(cAt);
      for (Integer dbs : sameC.keySet()){
        Hashtable<String,LipidParameterSet> sameDbs = sameC.get(dbs);
        if (!dbsRange.insideRange(dbs)){
          unusedSameC.put(dbs,sameDbs);
          continue;
        }
        Hashtable<String,Float> sameDbsDiffs = new Hashtable<String,Float>();
        for (String rt : sameDbs.keySet()){
          float[] input = new float[2];
          input[0] = cAt;
          input[1] = dbs;
          float difference = getRtFromRtString(rt)-optimizer.calculateFitValue(input);
          float tol = tolerance*optimizer.getMeanDeviation();
          if (COUNTER_SERIES_TOLERANCE>tol) tol = COUNTER_SERIES_TOLERANCE;
          if (Math.abs(difference)>tol) continue;
          ////if (Math.abs(difference)>(tolerance*optimizer.getMeanDeviation())) continue;
          sameDbsDiffs.put(rt, difference);
          if (difference<0) countNeg++;
          else countPos++;
        }
        if (sameDbsDiffs.size()>0) sameCDiffs.put(dbs, sameDbsDiffs);
      }
      if (sameCDiffs.size()>0) diffs.put(cAt, sameCDiffs);
      if (unusedSameC.size()>0) unusedNegs.put(cAt, unusedSameC);
    }
    // the should be applied only if 75% of the negatives are in one direction of the model
    if (countPos>=(3*countNeg) || countNeg>=(3*countPos) && (countPos+countNeg)>=MINIMUM_MEASUREMENTS){
      boolean usePos = true;
      if (countNeg>countPos) usePos = false;
      for (Integer cAt : diffs.keySet()){
        Hashtable<Integer,Hashtable<String,LipidParameterSet>> sameC = new Hashtable<Integer,Hashtable<String,LipidParameterSet>>();
        for (Integer dbs : diffs.get(cAt).keySet()){
          Hashtable<String,LipidParameterSet> sameDbs = new Hashtable<String,LipidParameterSet>();
          for (String rt : diffs.get(cAt).get(dbs).keySet()){
            float diff = diffs.get(cAt).get(dbs).get(rt);
            if ((diff>0 && usePos) || (diff<0 && !usePos)){
              sameDbs.put(rt, negatives.get(cAt).get(dbs).get(rt));
            }
          }
          if (sameDbs.size()>0) sameC.put(dbs, sameDbs);
        }
        if (sameC.size()>0) negs.put(cAt, sameC);
      }
    } else return null;
    results.add(negs);
    results.add(unusedNegs);
    return results;
  }
  
//  private boolean isCloserToCounterModel(float rt, float fitRt, float counterRt){
//    boolean closerToCounter = false;
//    if (!Float.isNaN(counterRt)){
//      if (Math.abs(rt-counterRt)<Math.abs(rt-fitRt)) closerToCounter = true;
//    }
//    return closerToCounter;
//  }
  
  /**
   * returns true if the retention time is closer to the counter model than the regular model
   * @param className name of lipid class
   * @param modName name of the modification
   * @param analyteName name of the analyte
   * @param rtString retention time as string
   * @param model the regular LM model
   * @param counterModel the counter LM model
   * @return true if the retention time is closer to the counter model than the regular model
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the library is not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there 
   * @throws LMException exception if model adaption is not possible
   */
  public static boolean fallsOnCounterModelSide(String className, String modName, String analyteName, String rtString, LevenbergMarquardtOptimizer model, LevenbergMarquardtOptimizer counterModel) throws RulesException, NoRuleException, IOException, SpectrummillParserException, LMException{
    String ruleName = StaticUtils.getRuleName(className, modName);
    Pattern cAtomsPattern =  Pattern.compile(RulesContainer.getCAtomsFromNamePattern(ruleName));
    Pattern dbsPattern =  Pattern.compile(RulesContainer.getDoubleBondsFromNamePattern(ruleName));
    Matcher cAtomsMatcher = cAtomsPattern.matcher(analyteName.split(LipidomicsConstants.CHAIN_MOD_SEPARATOR)[0]);
    if (!cAtomsMatcher.matches()) throw new RulesException("The analyte "+analyteName+" does not match the "+FragRuleParser.GENERAL_CATOMS_PARSE+" pattern \""+RulesContainer.getCAtomsFromNamePattern(ruleName)+"\" of the class "+ruleName+"!");
    Matcher dbsMatcher = dbsPattern.matcher(analyteName.split(LipidomicsConstants.CHAIN_MOD_SEPARATOR)[0]);
    int cAtoms = Integer.parseInt(cAtomsMatcher.group(1));
    if (!dbsMatcher.matches()) throw new RulesException("The analyte "+analyteName+" does not match the "+FragRuleParser.GENERAL_DBOND_PARSE+" pattern \""+RulesContainer.getDoubleBondsFromNamePattern(ruleName)+"\" of the class "+ruleName+"!");
    int dbs = Integer.parseInt(dbsMatcher.group(1));
    float rt = Float.valueOf(rtString);
    float[] cAndDbs = new float[2];
    cAndDbs[0] = cAtoms;
    cAndDbs[1] = dbs;
    float fitRt = model.calculateFitValue(cAndDbs);
    float counterRt = Float.NaN;
    if (counterModel!=null) counterRt = counterModel.calculateFitValue(cAndDbs);
    return fallsOnCounterModelSide(rt,fitRt,counterRt);
  }
  
  private static boolean fallsOnCounterModelSide(float rt, float fitRt, float counterRt){
    boolean onCounterSide = false;
    float percentToCount = 0.2f;
    if (!Float.isNaN(counterRt)){
      if (fitRt<counterRt){
        if (rt>counterRt || rt>(counterRt-(counterRt-fitRt)*percentToCount)) onCounterSide = true;
      } else {
        if (rt<counterRt || rt<(counterRt+(fitRt-counterRt)*percentToCount)) onCounterSide = true;
      }
    }
    return onCounterSide;
  }
  
  @SuppressWarnings("unchecked")
  private Vector<LevenbergMarquardtOptimizer> evaluateFilterAgainBasedOnCounterModel(Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> currentHits,
      Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> negatives,
      LevenbergMarquardtOptimizer model , LevenbergMarquardtOptimizer counter, boolean respectOh) throws LMException{
    int countNegatives = 0;
    @SuppressWarnings("rawtypes")
    Vector models = new Vector();
    Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> forModel = new Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>();
    for (Integer cAtoms: currentHits.keySet()){
      Hashtable<Integer,Hashtable<String,LipidParameterSet>> newSameC = new Hashtable<Integer,Hashtable<String,LipidParameterSet>>();
      Hashtable<Integer,Hashtable<String,LipidParameterSet>> negSameC = new Hashtable<Integer,Hashtable<String,LipidParameterSet>>();
      if (negatives.containsKey(cAtoms)) negSameC = negatives.get(cAtoms);
      for (Integer dbs : currentHits.get(cAtoms).keySet()){
        Hashtable<String,LipidParameterSet> newSameDbs = new Hashtable<String,LipidParameterSet>();
        Hashtable<String,LipidParameterSet> negSameDbs = new Hashtable<String,LipidParameterSet>();
        if (negSameC.containsKey(dbs)) negSameDbs = negSameC.get(dbs);
        for (String rtString : currentHits.get(cAtoms).get(dbs).keySet()){
          if (negatives.containsKey(cAtoms) && negatives.get(cAtoms).containsKey(dbs) && negatives.get(cAtoms).get(dbs).containsKey(rtString)){
            break;
          }
          float rt = getRtFromRtString(rtString);
          float[] cAndDbs = new float[2];
          cAndDbs[0] = cAtoms;
          cAndDbs[1] = dbs;
          float fitRt = model.calculateFitValue(cAndDbs);
          float counterRt = Float.NaN;
          if (counter!=null) counterRt = counter.calculateFitValue(cAndDbs);
          if (fallsOnCounterModelSide(rt,fitRt,counterRt)){
            negSameDbs.put(rtString, currentHits.get(cAtoms).get(dbs).get(rtString));
            countNegatives++;
          }else
            newSameDbs.put(rtString, currentHits.get(cAtoms).get(dbs).get(rtString));

        }
        if (newSameDbs.size()>0) newSameC.put(dbs, newSameDbs);
        if (negSameDbs.size()>0) negSameC.put(dbs, negSameDbs);
      }
      if (newSameC.size()>0) forModel.put(cAtoms, newSameC);
      if (negSameC.size()>0) negatives.put(cAtoms, negSameC);
    }
    if (countNegatives>0){
      LevenbergMarquardtOptimizer newModel = optimizeLMModel(forModel, null, respectOh);
      LevenbergMarquardtOptimizer newCounterModel = optimizeLMModel(negatives, null, respectOh);
      models.add(newModel);
      models.add(newCounterModel);
      models.add(negatives);
      return models;
    }else{
      models.add(model);
      models.add(counter);
      models.add(negatives);
      return models;
    }
  }
  
  private float getRtFromRtString(String rtString){
    String rt = new String (rtString);
    if (rt.indexOf("_")!=-1) rt = rt.substring(0,rt.indexOf("_"));
    return Float.valueOf(rt);
  }
  
  private float extractMinimumAcceptedDeviationValue(String className, Hashtable<String,Hashtable<String,LipidParameterSet>> data) throws RulesException, NoRuleException, IOException, SpectrummillParserException{
    Vector<Float> peakWidths = new Vector<Float>();
    for (Hashtable<String,LipidParameterSet> sameAnal : data.values()){
      for (LipidParameterSet set : sameAnal.values()){
        if (!set.isSuitableForRtProcessingHit(className)) continue;
        Vector<Vector<CgProbe>> probes = set.getIsotopicProbes();
        if (probes.size() == 0) continue;
        for (CgProbe probe: probes.get(0)){
          peakWidths.add((probe.UpperValley-probe.LowerValley)/60f);
        }
      }
    }
    float minDev = 0f;
    if (peakWidths.size()>0) minDev = Calculator.median(peakWidths)/(2f);
    return minDev;
  }
  
  
  /**
   * makes a selection between two equally matching adducts based on retention time
   * @param quantObjects the original quantitation instructions that were read from the Excel file
   * @return results after the selection was made
   * @throws RulesException specifies in detail which rule has been infringed
   * @throws NoRuleException thrown if the library is not there
   * @throws IOException general exception if there is something wrong about the file
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * @throws LMException exception if model adaption is not possible - should not be thrown, since only predicted models are used
   */
  public Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> chooseMoreLikelyOne(Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects) throws RulesException, NoRuleException, IOException, SpectrummillParserException, LMException {
    Set<String> affectedClasses = new HashSet<String>();
    Set<String> affectedMods = new HashSet<String>();
    Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> undecidedHits = new Hashtable<QuantVO,Hashtable<String,LipidParameterSet>>();
    Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> unprocessed = new Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>>();
    for (String aClass : this.results_.keySet()) {
      for (String anal : this.results_.get(aClass).keySet()) {
        for (String mod : this.results_.get(aClass).get(anal).keySet()) {
          for (String rt : this.results_.get(aClass).get(anal).get(mod).keySet()) {
            LipidParameterSet set = this.results_.get(aClass).get(anal).get(mod).get(rt);
            if (set.isChoseMoreLikelyRtWhenEqualMSn()) {
              affectedClasses.add(aClass);
              affectedMods.add(StaticUtils.getRuleName(aClass, set.getModificationName()));
              QuantVO quant = quantObjects.get(aClass).get(anal).get(mod);
              if (!undecidedHits.containsKey(quant))
                undecidedHits.put(quant, new Hashtable<String,LipidParameterSet>());
              undecidedHits.get(quant).put(rt, set);
            }
          }
        }
      }
    }
    for (String aClass : this.results_.keySet()) {
      if (affectedClasses.contains(aClass))
        unprocessed.put(aClass, this.results_.get(aClass));
    }
    correctByRetentionTimeSeries(unprocessed,ms2Removed_,adductInsensitiveRtFilter_);
    Vector<SharedMS1PeakVO> sharedPeaks = MSnPeakSeparator.detectSharedMS1PeakInstances(undecidedHits);
    Hashtable<String,LipidParameterSet> toRemove = new Hashtable<String,LipidParameterSet>();
    String idi;
    for (SharedMS1PeakVO shared : sharedPeaks) {
      if (shared.getPartners().size()<2 || shared.hasAnyPartnerDistinctFragments() || !shared.haveAllChooseOnRtSetToTrue())
        continue;
      boolean forAllModelsFound = true;
      float smallestDiff = Float.MAX_VALUE;
      String correct = null;
      for (SharedPeakContributionVO contr : shared.getPartners()) {
        LevenbergMarquardtOptimizer predictedModel = null;
        boolean respectOh = false;
        String ruleName = StaticUtils.getRuleName(contr.getQuantVO().getAnalyteClass(),contr.getQuantVO().getModName());
        if (adductInsensitiveRtFilter_.get(contr.getQuantVO().getAnalyteClass())) {
          if (predictedModels_.containsKey(contr.getQuantVO().getAnalyteClass())) {
            predictedModel = predictedModels_.get(contr.getQuantVO().getAnalyteClass());
            respectOh = respectOhs_.get(contr.getQuantVO().getAnalyteClass());
          }
        }else {
          if (predictedModels_.containsKey(ruleName)) {
            predictedModel = predictedModels_.get(ruleName);
            respectOh = respectOhs_.get(ruleName);
          }
        }
        if (predictedModel==null) {
          forAllModelsFound = false;
          break;
        }
        Pattern cAtomsPattern =  Pattern.compile(RulesContainer.getCAtomsFromNamePattern(ruleName));
        Pattern dbsPattern =  Pattern.compile(RulesContainer.getDoubleBondsFromNamePattern(ruleName));
        Matcher cAtomsMatcher = cAtomsPattern.matcher(contr.getQuantVO().getIdString().split(LipidomicsConstants.CHAIN_MOD_SEPARATOR)[0]);
        if (!cAtomsMatcher.matches()) throw new RulesException("The analyte "+contr.getQuantVO().getIdString()+" does not match the "+FragRuleParser.GENERAL_CATOMS_PARSE+" pattern \""+RulesContainer.getCAtomsFromNamePattern(ruleName)+"\" of the class "+ruleName+"!");
        int cAtoms = Integer.parseInt(cAtomsMatcher.group(1));
        Matcher dbsMatcher = dbsPattern.matcher(contr.getQuantVO().getIdString().split(LipidomicsConstants.CHAIN_MOD_SEPARATOR)[0]);
        if (!dbsMatcher.matches()) throw new RulesException("The analyte "+contr.getQuantVO().getIdString()+" does not match the "+FragRuleParser.GENERAL_DBOND_PARSE+" pattern \""+RulesContainer.getDoubleBondsFromNamePattern(ruleName)+"\" of the class "+ruleName+"!");
        int dbs = Integer.parseInt(dbsMatcher.group(1));
        float[] cAndDbs = new float[2];
        if (respectOh)
          cAndDbs = new float[3];
        cAndDbs[0] = cAtoms;
        cAndDbs[1] = dbs;
        if (respectOh)
          cAndDbs[2] = contr.getSet().getOhNumber();
        float fitRt = predictedModel.calculateFitValue(cAndDbs);
        float rt = Float.parseFloat(contr.getSet().getRt());
        float rtDiff = Math.abs(rt-fitRt);
        if (rtDiff<smallestDiff) {
          correct = OtherAdductChecker.getUniqueId(contr.getQuantVO().getAnalyteClass(), contr.getQuantVO().getIdString(), contr.getQuantVO().getModName(), contr.getSet().getRt());
        }
      }
      if (!forAllModelsFound)
        continue;
      for (SharedPeakContributionVO contr : shared.getPartners()) {
        idi = OtherAdductChecker.getUniqueId(contr.getQuantVO().getAnalyteClass(), contr.getQuantVO().getIdString(), contr.getQuantVO().getModName(), contr.getSet().getRt());
        if (!correct.equalsIgnoreCase(idi))
          toRemove.put(idi, contr.getSet());
        else
          contr.getSet().setPercentalSplit(100f);
      }
    }
    for (String id : toRemove.keySet()) {
      Vector<String> params = OtherAdductChecker.getParamsFromId(id);
      results_.get(params.get(0)).get(params.get(1)).get(params.get(2)).remove(params.get(3));
    }
    return results_;
  }
}
