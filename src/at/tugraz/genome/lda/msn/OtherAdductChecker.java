/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2019 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.Vector;

import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.exception.NoRuleException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.LipidomicsAnalyzer;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.QuantVO;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.quantification.CgException;
import at.tugraz.genome.maspectras.quantification.CgProbe;

/**
 * checks the results for other adducts and removes dependent results if they are not present
 * @author Juergen Hartler
 *
 */
public class OtherAdductChecker
{
  
  /** separator for an encoded string */
  private final static String SEP = "<->";
  
  
  /**
   * checks if the other required adducts are present to count this adduct as valid
   * @param results the detected results
   * @param unsplittedPeaks hash table containing unsplitted peaks; in case of a peak split, this hash is filled with unsplitted versions of the peak
   * @param quantObjects the original quantitation instructions that were read from the Excel file
   * @param analyzer object that holds MS data and can quantify fragments of interest
   * @param classSequence the MS levels of each lipid class
   * @return the results after removing adducts where the required other adducts are not present, and where the peaks were split correspondingly
   * @throws CgException exception that is thrown when there is something wrong with the fragment detection
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public static Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> checkTheResultsForOtherAdducts(
      Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> results,
      Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> unsplittedPeaks,
      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects, LipidomicsAnalyzer analyzer,
      LinkedHashMap<String,Integer> classSequence) throws CgException, LipidCombinameEncodingException{
    
    String ruleName;
    Hashtable<String,Boolean> requiresOtherAdducts = new Hashtable<String,Boolean>();
    Hashtable<String,Vector<String>> otherRequiredAdducts = new Hashtable<String,Vector<String>>();
    Hashtable<String,Boolean> allOtherAdductsRequired = new Hashtable<String,Boolean>();
    Hashtable<String,Float> timeTolerance = new Hashtable<String,Float>();
    boolean requiresOther;
    boolean forceAdductValidity = false;
    String quantId;
    Hashtable<String,LipidParameterSet> correct = new Hashtable<String,LipidParameterSet>();
    Hashtable<String,LipidParameterSet> toRemove = new Hashtable<String,LipidParameterSet>();
    Hashtable<String,String> interestingQuantVOs = new Hashtable<String,String>();
    Hashtable<String,Boolean> affectedMods = new Hashtable<String,Boolean>();

    //build the other adduct requirement hash
    for (String lipidClass : results.keySet()) {
      for (String analyte : results.get(lipidClass).keySet()) {
        for (String mod : results.get(lipidClass).get(analyte).keySet()) {
          requiresOther = false;
          ruleName = StaticUtils.getRuleName(lipidClass,mod);
          try{
            if (RulesContainer.requiresOtherValidAdduct(ruleName)) {
              requiresOther = true;
              requiresOtherAdducts.put(ruleName, requiresOther);
              otherRequiredAdducts.put(ruleName,RulesContainer.getOtherRequiredAdducts(ruleName));
              allOtherAdductsRequired.put(ruleName,RulesContainer.areAllOtherAdductsRequired(ruleName));
              timeTolerance.put(ruleName,RulesContainer.getOtherTimeTolerance(ruleName));
              forceAdductValidity = RulesContainer.forceOtherAdductValidity(ruleName);
              affectedMods.put(ruleName, forceAdductValidity);
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
    
    for (String lipidClass : results.keySet()) {
      for (String analyte : results.get(lipidClass).keySet()) {
        for (String mod : results.get(lipidClass).get(analyte).keySet()) {
          ruleName = StaticUtils.getRuleName(lipidClass,mod);
          Hashtable<String,LipidParameterSet> sameMod = results.get(lipidClass).get(analyte).get(mod);
          if (!requiresOtherAdducts.containsKey(ruleName) || requiresOtherAdducts.get(ruleName)==false)
            continue;
          for (String rt : sameMod.keySet()) {
            boolean oneFound = false;
            boolean allFound = true;
            for (String otherMod : otherRequiredAdducts.get(ruleName)) {
              //System.out.println(lipidClass+analyte+"_"+mod+"_"+rt+":  "+otherMod+"-"+isOtherModPresentInResults(lipidClass,analyte,otherMod,rt,timeTolerance,results));
              if (isOtherModPresentInResults(lipidClass,analyte,otherMod,rt,timeTolerance.get(ruleName),results))
                oneFound = true;
              else
                allFound = false;
            }
            quantId = getUniqueId(lipidClass,analyte,mod,"");
            interestingQuantVOs.put(quantId, quantId);
            // the other adducts were found
            if (oneFound && (!allOtherAdductsRequired.get(ruleName) || allFound)) 
              correct.put(getUniqueId(lipidClass,analyte,mod,rt),sameMod.get(rt));
            // the other adducts were not found
            else
              toRemove.put(getUniqueId(lipidClass,analyte,mod,rt),sameMod.get(rt));
          }
        }
      }
    }
    
    //now enforce the adduct or split peaks having different adducts
    //TODO: these lines have to be tested on a bigger data set!!!
    Vector<QuantVO> allQuants;
    for (String lClass : quantObjects.keySet()) {
      for (String analyte : quantObjects.get(lClass).keySet()){
        for (String mod : quantObjects.get(lClass).get(analyte).keySet()) {
          QuantVO quantSet = quantObjects.get(lClass).get(analyte).get(mod);
          allQuants = new Vector<QuantVO>();
          boolean isAffected = false;
          allQuants.add(quantSet);
          allQuants.addAll(quantSet.getOtherIsobaricSpecies());
          //check whether this quantVO is shared with something requiring other adducts, and if so, if the other one was found
          for (QuantVO one : allQuants) {
            if (!affectedMods.containsKey(StaticUtils.getRuleName(one.getAnalyteClass(),one.getModName())))
              continue;
            if (interestingQuantVOs.containsKey(getUniqueId(lClass,analyte,mod,"")))
              isAffected = true;
          }
          if (!isAffected)
            continue;
          Hashtable<String,Hashtable<String,LipidParameterSet>> sameRt = getGroupsOfSameRt(allQuants,results);
          String uniqueId;
          for (String rt : sameRt.keySet()) {
            Hashtable<String,LipidParameterSet> same = sameRt.get(rt);
            //nothing has to be done where there is no uncertainty
            if (same.size()<2)
              continue;
            Hashtable<String,LipidParameterSet> force = new Hashtable<String,LipidParameterSet>();
            boolean ms2InformationPresent = false;
            for (String rule : same.keySet()) {
              LipidParameterSet set = same.get(rule);
              String lipClass = rule.substring(0,rule.indexOf("_"));
              if (set instanceof LipidomicsMSnSet)
                ms2InformationPresent = true;
              if (!affectedMods.containsKey(rule) || !affectedMods.get(rule))
                continue;
              uniqueId = getUniqueId(lipClass,set.getNameStringWithoutRt(),set.getModificationName(),rt);
              if (correct.containsKey(uniqueId)) {
                force.put(uniqueId, set);
                //System.out.println("I have to force the following ones: "+lClass+analyte+"_"+mod+"_"+rt);
              }
            }
            // if there are hits that have to be forced - all the other ones have to be removed
            Hashtable<String,LipidParameterSet> forSplit = new Hashtable<String,LipidParameterSet>();
            if (force.size()>0) {
              forSplit = force;
              for (String rule : same.keySet()) {
                String lipClass = rule.substring(0,rule.indexOf("_"));
                LipidParameterSet set = same.get(rule);
                uniqueId = getUniqueId(lipClass,set.getNameStringWithoutRt(),set.getModificationName(),rt);
                if (!forSplit.containsKey(uniqueId))
                  toRemove.put(uniqueId, set);
              }
            // if there are no hits to force, all other ones have to be added to the split that are correct
            } else {
              for (String rule : same.keySet()) {
                String lipClass = rule.substring(0,rule.indexOf("_"));
                LipidParameterSet set = same.get(rule);
                uniqueId = getUniqueId(lipClass,set.getNameStringWithoutRt(),set.getModificationName(),rt);
                if (!toRemove.containsKey(uniqueId))
                  forSplit.put(uniqueId, set);
              }
            }
            
            if (!ms2InformationPresent)
              continue;
            if (forSplit.size()<2)
              continue;
            Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> hitsWithQuant = new Hashtable<QuantVO,Hashtable<String,LipidParameterSet>>();
            Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> peaksBeforeSplit = new Hashtable<QuantVO,Hashtable<String,LipidParameterSet>>();
            String lipClass = null;
            for (String id : forSplit.keySet()) {
              LipidParameterSet set = forSplit.get(id);
              lipClass = getLipidClassFromUniqueId(id);
              QuantVO quant = null;
              for (QuantVO one : allQuants) {
                if (lipClass.equalsIgnoreCase(one.getAnalyteClass()) && set.getNameStringWithoutRt().equalsIgnoreCase(one.getIdString()) &&
                  set.getModificationName().equalsIgnoreCase(one.getModName())) {
                  quant = one;
                  break;
                }
              }
              if (!hitsWithQuant.containsKey(quant))
                hitsWithQuant.put(quant, new Hashtable<String,LipidParameterSet>());
              hitsWithQuant.get(quant).put(set.getRt(), set);
              peaksBeforeSplit.put(quant, new Hashtable<String,LipidParameterSet>());
            }
            MSnPeakSeparator separator = new MSnPeakSeparator(hitsWithQuant, peaksBeforeSplit, analyzer,  classSequence.get(lipClass),new HashSet<String>());
            Hashtable<QuantVO,Hashtable<String,LipidParameterSet>> hitsAccordingToQuant = separator.disentagleSharedMS1Peaks();
            //if the hit was removed in the splitting process, add it for removal here
            for (QuantVO quant : hitsWithQuant.keySet()) {
              for (String ret : hitsWithQuant.get(quant).keySet()) {
                if (!hitsAccordingToQuant.containsKey(quant) || !hitsAccordingToQuant.get(quant).containsKey(ret)) {
                  toRemove.put(getUniqueId(quant.getAnalyteClass(),quant.getIdString(),quant.getModName(),ret), hitsWithQuant.get(quant).get(ret));
                }
              }
            }
            //replace the unsplitted hits with the splitted ones
            for (QuantVO quant : hitsAccordingToQuant.keySet()) {
              for (String ret : hitsAccordingToQuant.get(quant).keySet()) {
                LipidParameterSet set = hitsAccordingToQuant.get(quant).get(ret);
                results.get(quant.getAnalyteClass()).get(quant.getIdString()).get(set.getModificationName()).put(ret, set);
              }
            }
            for (QuantVO quant : peaksBeforeSplit.keySet()) {
              for (String ret : peaksBeforeSplit.get(quant).keySet()) {
                LipidParameterSet set = peaksBeforeSplit.get(quant).get(ret);
                Hashtable<String,LipidParameterSet> unsplitted = new Hashtable<String,LipidParameterSet>();
                if (unsplittedPeaks.get(quant.getAnalyteClass()).get(quant.getIdString()).containsKey(set.getModificationName()))
                  unsplitted = unsplittedPeaks.get(quant.getAnalyteClass()).get(quant.getIdString()).get(set.getModificationName());
                unsplitted.put(ret, set);
                unsplittedPeaks.get(quant.getAnalyteClass()).get(quant.getIdString()).put(set.getModificationName(), unsplitted);
              }
            }              
          }
        }
      }
    }
    
    //finally, remove all the ones where no other mandatory adducts were found
    for (String id : toRemove.keySet()) {
      Vector<String> params = getParamsFromId(id);
      results.get(params.get(0)).get(params.get(1)).get(params.get(2)).remove(params.get(3));
    }
    return results;
  }
  
  
  /**
   * searches the results whether another required modification is present in the results
   * @param lClass the lipid class to look for
   * @param analyte the name of the analyte
   * @param mod the name of the modification
   * @param rtString the string representation of the retention time to look for
   * @param tolerance the allowed retention time tolerance
   * @param results the detected results, i.e. the place to look for the other adducts
   * @return true when the other required modifications were found
   */
  private static boolean isOtherModPresentInResults(String lClass, String analyte, String mod, String rtString, float tolerance,
      Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> results) {
    if (!results.containsKey(lClass) || !results.get(lClass).containsKey(analyte) || !results.get(lClass).get(analyte).containsKey(mod))
      return false;
    boolean found = false;
    float rt = Float.parseFloat(rtString);
    float startRt = rt-tolerance;
    float stopRt = rt+tolerance;
    float foundRt;
    //iterate over the found RTs and check which ones are within the tolerance
    for (String foundRtString : results.get(lClass).get(analyte).get(mod).keySet()) {
      foundRt = Float.parseFloat(foundRtString);
      if (startRt<=foundRt && foundRt<=stopRt) {
        found = true;
        break;
      }
    }    
    return found;
  }
 
  
  /**
   * creates a unique id as key for hash tables consisting of the method parameters
   * @param lClass lipid class
   * @param anal name of the analyte
   * @param mod name of the modification
   * @param rt string representation of the retention time
   * @return unique id as key for hash tables consisting of the method parameters
   */
  private static String getUniqueId(String lClass, String anal, String mod, String rt) {
    return lClass+SEP+anal+SEP+mod+SEP+rt;
  }
  
  
  /**
   * parses the lipid class out of the unique id -> see getUniqueId
   * @param uniqueId the unique id string
   * @return the name of the lipid class
   */
  private static String getLipidClassFromUniqueId(String uniqueId) {
    return uniqueId.substring(0, uniqueId.indexOf(SEP));
  }
  
  
  /**
   * checks for isobaric/isomeric detections that originate from the same peak (same retention time)
   * @param allQuants the isobaric/isomeric quantiation instructions
   * @param results the detected results
   * @return isobaric/isomeric detections that originate from the same peak (same retention time) - first key: retention time; second key: rule name consisting of lipid class and modification
   */
  private static Hashtable<String,Hashtable<String,LipidParameterSet>> getGroupsOfSameRt(Vector<QuantVO> allQuants,
      Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> results){
    Hashtable<String,Hashtable<String,LipidParameterSet>> sameRts = new Hashtable<String,Hashtable<String,LipidParameterSet>>();
    for (QuantVO quant : allQuants) {
      if (!results.containsKey(quant.getAnalyteClass()) || !results.get(quant.getAnalyteClass()).containsKey(quant.getIdString())||
          !results.get(quant.getAnalyteClass()).get(quant.getIdString()).containsKey(quant.getModName()))
      continue;
      for (String rt : results.get(quant.getAnalyteClass()).get(quant.getIdString()).get(quant.getModName()).keySet()) {
        LipidParameterSet set = results.get(quant.getAnalyteClass()).get(quant.getIdString()).get(quant.getModName()).get(rt);
        if (!sameRts.containsKey(rt))
          sameRts.put(rt, new Hashtable<String,LipidParameterSet>());
        sameRts.get(rt).put(StaticUtils.getRuleName(quant.getAnalyteClass(),quant.getModName()),set);          
      }
    }
    return sameRts;
  }
  
  /**
   * returns the individual parameters of a unique id that was used for the hash tables - see getUniqueId
   * @param id the unique id that was used for the hash tables
   * @return the individual parameters get(0): lipid class; get(1): analyte name; get(2): modification name; get(3): retention time
   */
  private static Vector<String> getParamsFromId(String id){
    Vector<String> params = new Vector<String>();
    String rest = new String(id);
    while (rest.length()>0) {
      if (rest.indexOf(SEP)!=-1) {
        params.add(rest.substring(0,rest.indexOf(SEP)));
        rest = rest.substring(rest.indexOf(SEP)+SEP.length());
      } else {
        params.add(rest);
        rest = "";
      }
    }
    return params;
  }
  
  
  /**
   * removes peaks that fall below the base peak cutoff
   * @param results the results where lower entities shall be removed
   * @param threshold the threshold for the removal
   */
  public static void removePeaksThatFallBelowTheBasepeakCutoff(Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> results , float threshold) {
    Hashtable<String,LipidParameterSet> toRemove = new Hashtable<String,LipidParameterSet>();
    for (String anClass : results.keySet()) {
      for (String analyte : results.get(anClass).keySet()) {
        for (String mod : results.get(anClass).get(analyte).keySet()) {
          for (String rt : results.get(anClass).get(analyte).get(mod).keySet()) {
            LipidParameterSet set = results.get(anClass).get(analyte).get(mod).get(rt);
            boolean oneOk = false;
            for (CgProbe probe : set.getIsotopicProbes().get(0)){
              if (probe.Area>threshold)
                oneOk = true;
            }
            if (!oneOk)
              toRemove.put(getUniqueId(anClass,analyte,mod,rt),set);
          }
        }
      }
    }
    for (String id : toRemove.keySet()) {
      Vector<String> params = getParamsFromId(id);
      results.get(params.get(0)).get(params.get(1)).get(params.get(2)).remove(params.get(3));
    }
  }
  
}
