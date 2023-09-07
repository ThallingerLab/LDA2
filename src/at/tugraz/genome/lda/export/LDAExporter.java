/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2018 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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

package at.tugraz.genome.lda.export;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.ExportException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.exception.RetentionTimeGroupingException;
import at.tugraz.genome.lda.export.vos.EvidenceBase;
import at.tugraz.genome.lda.export.vos.EvidenceVO;
import at.tugraz.genome.lda.export.vos.FeatureVO;
import at.tugraz.genome.lda.export.vos.SpeciesExportVO;
import at.tugraz.genome.lda.export.vos.SummaryVO;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.msn.vos.IntensityPositionVO;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.DoubleBondPositionVO;
import at.tugraz.genome.lda.vos.DoubleStringVO;
import at.tugraz.genome.lda.vos.ResultAreaVO;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.voutils.GeneralComparator;

/**
 * This class is the base for exporting information based on selections in the heat maps plus extracting detailed
 * information stored in the referenced LDA result files
 * 
 * @author Juergen Hartler
 *
 */

public abstract class LDAExporter
{
  
  /**
   * extracts original information that belongs exactly to one data point (rectangle) in the heat map
   * @param all all detected results of a certain experiment and lipid class
   * @param areaVO the value object corresponding to one data point (rectangle) in the heat map
   * @param msnOnly if true, only species verified by MSn evidence are exported
   * @return the original information relevant to this one data point (rectangle)
   */
  protected static Hashtable<String,Vector<LipidParameterSet>> getRelevantOriginalResults(Vector<LipidParameterSet> all, ResultAreaVO areaVO, boolean msnOnly){
    Hashtable<String,Vector<LipidParameterSet>> results = new Hashtable<String,Vector<LipidParameterSet>>();
    for (LipidParameterSet set : all){
      if (msnOnly && !(set instanceof LipidomicsMSnSet))
        continue;
      //the next line is for the export of MSn data only
//      if (!(set instanceof LipidomicsMSnSet))
//        continue;
      if (!areaVO.getMoleculeNameWoRT().equalsIgnoreCase(set.getNameStringWithoutRt()))
        continue;
      if (!areaVO.belongsRtToThisAreaVO(set.getRt(), set.getModificationName()))
        continue;
      Vector<LipidParameterSet> ofOneMod = new Vector<LipidParameterSet>();
      if (results.containsKey(set.getModificationName()))
        ofOneMod = results.get(set.getModificationName());
      ofOneMod.add(set);
      results.put(set.getModificationName(), ofOneMod);
    }
    return results;
  }
  
  
  /**
   * checks whether this species is in vector (respects potential permutations)
   * @param species the species to check
   * @param molecularSpecies the vector containing other species
   * @return the String representation of the species in the list
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  protected static String isSpeciesAlreadyInList(String species, Vector<String> molecularSpecies) throws LipidCombinameEncodingException{
    for (String other : molecularSpecies){
      if (species.equalsIgnoreCase(other) || StaticUtils.isAPermutedVersion(species, other,LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)){
        return other;
      } 
    }  
    return null;
  }

  /**
   * creates a sorted hash map where the key is the representation with position information (if
   * at least found in more than half of the MS-detections and no contradicting evidence exists),
   * and the value is the used representation without position information
   * @param molecularSpecies the representations without position
   * @param molSpeciesCount the total count of detections of this species without position
   * @param positionCount a hash table containing the counts of this species; first key: without
   * position representation, second key: representation with position information; value: number
   * of occurrences of this position representation
   * @param ambiguous a hash table containing the counts where more than one possible position 
   * assignment was made; first key: without position representation, second key: representation
   * with position information; value: number of occurrences of this position representation
   * @return sorted hash map where the key is the representation with position information and the
   * value is the used representation without position information
   */
  protected static HashMap<String,String> replaceSpeciesByPositionAccordingToEvidence(Vector<String> molecularSpecies,
      Hashtable<String,Integer> molSpeciesCount, Hashtable<String,Hashtable<String,Integer>> positionCount,
      Hashtable<String,Vector<String>> ambiguous){
    HashMap<String,String> newValues = new HashMap<String,String>();
    for (String woPosition : molecularSpecies){
      int totalCount = molSpeciesCount.get(woPosition);
      Hashtable<String,Integer> positionsFound = positionCount.get(woPosition);
      boolean usePositionIdentification = false;
      //use position if there is no contradicting evidence
      if (positionsFound.size()==1){
        ////the second check was removed for sphingolipid export
        int detections = positionsFound.values().iterator().next();
        //this is for ambiguous identifications
        if (ambiguous!=null && ambiguous.containsKey(woPosition)){
          boolean allDetected = true;
          String detected = positionsFound.keySet().iterator().next();
          for (String ambi : ambiguous.get(woPosition)){
            boolean oneTimeCorrect = false;
            String[] possibilities = ambi.split(";");
            for (String onePossibility : possibilities){
              if (onePossibility.equalsIgnoreCase(detected)){
                oneTimeCorrect = true;
                break;
              }
            }
            if (oneTimeCorrect)
              detections++;
            else{
              allDetected = false;
              break;
            }
          }
          if (!allDetected){
            detections = 0;
          }
        }
        //use a position when it was found more than once, and at least more than half of the times of the
        //detected cases
        ////TODO: the next line is only excluded for JOE - I am not sure whether I actually need this removal!!!!!
//        if (detections>1 && detections>(totalCount/2))
        //the line before, is the last line for the removal, for sphingolipid export  
          usePositionIdentification = true;
      }
        
      if (usePositionIdentification){
        newValues.put(positionsFound.keySet().iterator().next(),woPosition);
      }else{
        newValues.put(woPosition,woPosition);
      }
      
    }
    return newValues;
  }
  

  
  
  /**
   * extracts the combined information of several experiments according to selections in the heat maps plus detailed data from LDA result files
   * @param speciesType structural level of data export (lipid species, chain level, position level - for details see LipidomicsConstants.EXPORT_ANALYTE_TYPE)
   * @param exportDoubleBondPositionsForClass true when double bond positions shall be exported for the lipid class of this analyte
   * @param extractFeatures should the features information be included in the exported data
   * @param currentSummaryId a running unique identifier
   * @param currentFeatureId a running unique identifier
   * @param extractEvidence should the evidence information information be included in the exported data (only possible when features are included)
   * @param currentEvidenceId a running unique identifier
   * @param currentEvGroupingId a running identifier for evidence originating from the same spectra
   * @param isRtGrouped is there only one entry for each species, or are the species grouped by retention time
   * @param adductsSorted key: the adducts sorted in consecutive manner starting with the strongest representative; value: contains this adduct position information 
   * @param expNames the sorted experiment names
   * @param expsOfGroup key: group name; value: experiments belonging to this group
   * @param molName the name of the molecule (inclusively an RT grouping parameter)
   * @param resultsMol the values according to the heat map selection
   * @param relevantOriginals the relevant results from the LDA Excel file - when these values are null, the speciesType must be LipidomicsConstants.EXPORT_ANALYTE_TYPE_SPECIES, and extractFeatures, extractEvidence and exportDoubleBondPositionsForClass must be false 
   * @param maxIsotope the maximum isotope that will be used for the export
   * @param faHydroxyEncoding the OH encodings of the FA moiety
   * @param lcbHydroxyEncoding the OH encodings of the LCB moiety
   * @return value object containing all exportable information
   * @throws ExportException when there is something wrong
   * @throws SpectrummillParserException when there are elements missing in the elementconfig.xml
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  @SuppressWarnings("unchecked")
  protected static SpeciesExportVO extractExportableSummaryInformation(short speciesType, boolean exportDoubleBondPositionsForClass, boolean extractFeatures, int currentSummaryId,
      int currentFeatureId, boolean extractEvidence, int currentEvidenceId, int currentEvGroupingId, boolean isRtGrouped,
      LinkedHashMap<String,Boolean> adductsSorted, Vector<String> expNames, LinkedHashMap<String,Vector<String>> expsOfGroup,
      String molName, Hashtable<String,Vector<Double>> resultsMol, Hashtable<String,Hashtable<String,Vector<LipidParameterSet>>> relevantOriginals,
      int maxIsotope, HydroxyEncoding faHydroxyEncoding, HydroxyEncoding lcbHydroxyEncoding)
          throws ExportException, SpectrummillParserException, LipidCombinameEncodingException,RetentionTimeGroupingException{
    if (extractEvidence && !extractFeatures)
      throw new ExportException("It is not possible to extract the evidence without extracting the features");
    String chemFormula = "";
    Double neutralMassTheoretical = null;
    Vector<DoubleBondPositionVO> allAssignedDoubleBondPositions = new Vector<DoubleBondPositionVO>();
    //reading information that is equal for all the identifications of one analyte
    if (relevantOriginals!=null){
      for (Hashtable<String,Vector<LipidParameterSet>> origs : relevantOriginals.values()){
        if (origs.size()==0)
          continue;
        if (chemFormula.length()>0)
          break;
        for (Vector<LipidParameterSet> sets : origs.values()){
          if (sets.size()==0)
            continue;
          if (chemFormula.length()>0)
            break;
          for (LipidParameterSet set : sets){
            try {
              Hashtable<String,Integer> categorized = StaticUtils.categorizeFormula(set.getAnalyteFormula());
              chemFormula = StaticUtils.getFormulaInHillNotation(categorized, false);
              neutralMassTheoretical = Settings.getElementParser().calculateTheoreticalMass(StaticUtils.getFormulaInHillNotation(categorized, true), false);
              break;
            }catch (ChemicalFormulaException e) {e.printStackTrace();}
          }
        }
      }
    }
    Vector<SummaryVO> speciesSummaries = new Vector<SummaryVO>();
    String molNameWORt = new String(molName);
    if (isRtGrouped){
      molNameWORt = molName.substring(0,molName.lastIndexOf("_"));
    }
    
    
    //this generates the feature objects
    LinkedHashMap<String,FeatureVO> features = null;
    if (extractFeatures){
      features = new LinkedHashMap<String,FeatureVO>();
      for (String mod : adductsSorted.keySet()){
        boolean foundOneHit = false;
        double totalArea = 0d;
        double areaTimesMz = 0d;
        double weightedMz = 0d;
        double areaTimesRt = 0d;
        double weightedRt = 0d;
        int charge = 0;
        float startRt = Float.MAX_VALUE;
        float stopRt = 0f;
        Vector<Double> areas = new Vector<Double>();
        for (String expName : expNames){
          if (!relevantOriginals.containsKey(expName) || !relevantOriginals.get(expName).containsKey(mod)){
            areas.add(null);
            continue;
          }
          Vector<LipidParameterSet> sets = relevantOriginals.get(expName).get(mod);
          double areaOfAssay = 0d;
          for (LipidParameterSet set : sets){
            if (charge==0)
              charge = set.getCharge();
            double areaOfOnePeak = set.getArea(maxIsotope);
            //calculate for each set a weighted mean
            double areaOfSet = 0d;
            double areaTimesMzOfSet = 0d;
            double areaTimesRtOfSet = 0d;
            boolean lowerRtHardLimitSet = false;
            boolean upperRtHardLimitSet = false;
            if (set.getLowerRtHardLimit()>0f){
              if (set.getLowerRtHardLimit()<startRt)
                startRt = set.getLowerRtHardLimit();
              lowerRtHardLimitSet = true;
            }
            if (set.getUpperRtHardLimit()>0f){
              if (set.getUpperRtHardLimit()>stopRt)
                stopRt = set.getUpperRtHardLimit();
              upperRtHardLimitSet = true;
            }
            for (CgProbe probe : set.getIsotopicProbes().get(0)){
              areaOfSet += (double)probe.Area;
              areaTimesMzOfSet += (double)probe.Area*(double)probe.Mz;
              areaTimesRtOfSet += (double)probe.Area*(double)probe.Peak;
              if (!lowerRtHardLimitSet && probe.LowerValley<startRt)
                startRt = probe.LowerValley;
              if (!upperRtHardLimitSet && probe.UpperValley>stopRt)
                stopRt = probe.UpperValley;
            }
            double weightedMeanArea = areaTimesMzOfSet/areaOfSet;
            double weightedMeanRt = areaTimesRtOfSet/areaOfSet;
            totalArea += areaOfOnePeak;
            areaOfAssay += areaOfOnePeak;
            areaTimesMz += weightedMeanArea*areaOfOnePeak;
            areaTimesRt += weightedMeanRt*areaOfOnePeak;
          }
          if (sets.size()>0){
            foundOneHit = true;
            weightedMz = areaTimesMz/totalArea;
            weightedRt = areaTimesRt/totalArea;
          }
          areas.add(areaOfAssay>0d ? areaOfAssay : null);
        }
        if (foundOneHit){
          FeatureVO feature = new FeatureVO(currentFeatureId,mod,weightedMz,charge,weightedRt,
              (double)startRt,(double)stopRt,areas);
          currentFeatureId++;
          features.put(mod, feature);
        }
      }
    }

    if (speciesType==LipidomicsConstants.EXPORT_ANALYTE_TYPE_SPECIES){
      Hashtable<String,Double> areas = calculateRelativeAreas(null,resultsMol,null);
      String strongestExp = getStrongestExp(areas);
      //detect the found modifications for this hit
      //detect the reliability of the hit 2 for MS/MS verification; 3 for MS1 only
      //extract the retention times of the highest peaks of every modification; rtsOfMods: first key modification; second key: experiment
      int mzTabReliability = 3;
      float rTime = 0f;
      Vector<String> adducts = null;
      Hashtable<String,Hashtable<String,Double>> rtsOfMods = new Hashtable<String,Hashtable<String,Double>>();
      //the reliability of the evidence of every modification; evidenceReliabilityOfMods: first key: species second key modification; third key: experiment
      Hashtable<String,Hashtable<String,Hashtable<String,Short>>> evidenceReliabilityOfMods = new Hashtable<String,Hashtable<String,Hashtable<String,Short>>>();
      if (relevantOriginals!=null){
        Hashtable<String,String> foundMods = new Hashtable<String,String>();
        short evidence;
        short otherEvidence;
        for (String exp : relevantOriginals.keySet()){
          Hashtable<String,Vector<LipidParameterSet>> modParams = relevantOriginals.get(exp);
          
          if (exportDoubleBondPositionsForClass) {
            for (Vector<LipidParameterSet> sets : modParams.values()) {
              for (LipidParameterSet set : sets) {
                if (set.hasOmegaInformation()) {
                  allAssignedDoubleBondPositions.addAll(StaticUtils.getAssignedDoubleBondPositions(set.getOmegaInformation()));
                }
              }
            }
          }
          
          for (String mod : modParams.keySet()){
            foundMods.put(mod, mod);
            evidence = SummaryVO.EVIDENCE_MS1_ONLY;
            for (LipidParameterSet set : modParams.get(mod)){
              if (mzTabReliability == 3){
                if (set instanceof LipidomicsMSnSet && ((LipidomicsMSnSet)set).getStatus()>=LipidomicsMSnSet.HEAD_GROUP_DETECTED)
                  mzTabReliability = 2;
              }
              otherEvidence = StaticUtils.determineEvidenceStateOfHit(set);
              if (otherEvidence>evidence)
                evidence = otherEvidence;
            }
            storeEvidenceToHash(evidence, evidenceReliabilityOfMods, molName, mod, exp);
            @SuppressWarnings("rawtypes")
            Vector rtAndArea = getRtOfHighestZeroIsoArea(modParams.get(mod));
            if (!rtsOfMods.containsKey(mod))
              rtsOfMods.put(mod, new Hashtable<String,Double>());
            rtsOfMods.get(mod).put(exp, ((Double)rtAndArea.get(0)/60d));
          }
        }
        //detect the Rt of the strongest identification
        if (strongestExp!=null){
          Vector<LipidParameterSet> allParams = new Vector<LipidParameterSet>();
          for (Vector<LipidParameterSet> params : relevantOriginals.get(strongestExp).values()){
            allParams.addAll(params);
          }
          @SuppressWarnings("rawtypes")
          Vector rtAndArea = getRtOfHighestZeroIsoArea(allParams);
          rTime = ((Double)rtAndArea.get(0)).floatValue();
        }
        adducts = getSortedModifications(adductsSorted,foundMods);
      }
      Vector<Integer> featureRefs = getFeatureRefs(molName,adducts,features);
      String displayString = null;
      if (exportDoubleBondPositionsForClass) {
        for (DoubleBondPositionVO assignedDoubleBondPosition : allAssignedDoubleBondPositions) {
          String assignedSpecies = assignedDoubleBondPosition.getDoubleBondPositionsHumanReadable(speciesType);
          if (displayString != null && !displayString.equals(assignedSpecies)) {
            System.out.println(String.format("Analytes with differing assigned double bond positions have been grouped together (%s and %s)."
                + "Consider lowering the setting for 'Show hits with different RT separately' in the 'Selection' tab.", displayString, assignedSpecies));
            displayString = null;
            break;    
//            throw new RetentionTimeGroupingException(
//                "<html>Analytes with differing assigned double bond positions have been falsely grouped together.<br/>"
//                + "Consider lowering the setting for 'Show hits with different RT separately' in the 'Selection' tab.</html>");
          }
          displayString = assignedSpecies;
        }
      }
      SummaryVO sumVO = new SummaryVO(currentSummaryId,molName,displayString,featureRefs,chemFormula,neutralMassTheoretical,rTime,
          adducts,mzTabReliability,areas,rtsOfMods,evidenceReliabilityOfMods.get(molName),expsOfGroup);
      currentSummaryId++;
      speciesSummaries.add(sumVO);
    } else if (speciesType==LipidomicsConstants.EXPORT_ANALYTE_TYPE_CHAIN || speciesType==LipidomicsConstants.EXPORT_ANALYTE_TYPE_POSITION){
      //this is a Vector containing species names only if there is no molecular species found in this result file
      Vector<String> species = new Vector<String>();
      int mzTabReliabilityOfSumSpecies = 3;
      //this vector contains the detected molecular species without position detection
      Vector<String> molecularSpecies = new Vector<String>();
      //this hashtable is counting how often a certain position was found - the first key is the the molecular species without position assignment,
      //the second key is the molecular species with a position assigned, and the value is the number of detections
      Hashtable<String,Hashtable<String,Integer>> positionCount = new Hashtable<String,Hashtable<String,Integer>>();
      //this hashtable is counting how often a molecular species was found
      Hashtable<String,Integer> molSpeciesCount = new Hashtable<String,Integer>();
      //this hashtable contains ambiguous identifications
      Hashtable<String,Vector<String>> ambiguousOnes = new Hashtable<String,Vector<String>>();
      HashMap<String,String> positionToWithoutHash = null; 
      
      //first key: experiment name; second key species name; the value is the percental split;
      Hashtable<String,Hashtable<String,Double>> percentalSplits = new Hashtable<String,Hashtable<String,Double>>();
      Hashtable<String,Hashtable<String,Float>> highestRts = new Hashtable<String,Hashtable<String,Float>>(); 
      //first key: species; second key and value is the modification
      Hashtable<String,Hashtable<String,String>> modsOfSpecies = new Hashtable<String,Hashtable<String,String>>();
      
      //the retention times of the highest peaks of every modification; rtsOfMods: first key: species second key modification; third key: experiment
      Hashtable<String,Hashtable<String,Hashtable<String,Double>>> rtsOfMods = new Hashtable<String,Hashtable<String,Hashtable<String,Double>>>();

      //the reliability of the evidence of every modification; evidenceReliabilityOfMods: first key: species second key modification; third key: experiment
      Hashtable<String,Hashtable<String,Hashtable<String,Short>>> evidenceReliabilityOfMods = new Hashtable<String,Hashtable<String,Hashtable<String,Short>>>();
      
      //this is for the summary information
      for (String expName : expNames){
        if (!relevantOriginals.containsKey(expName))
          continue;
        Hashtable<String,Vector<LipidParameterSet>> allMods = relevantOriginals.get(expName);
        Vector<LipidParameterSet> allSets = new Vector<LipidParameterSet>();
        for ( Vector<LipidParameterSet> sets : allMods.values()) {
          allSets.addAll(sets);
        } 
        if (exportDoubleBondPositionsForClass) {
          for (LipidParameterSet set : allSets) {
            if (set.hasOmegaInformation()) {
              allAssignedDoubleBondPositions.addAll(StaticUtils.getAssignedDoubleBondPositions(set.getOmegaInformation()));
            }
          }
        }
        Hashtable<String,Float> highestPeakArea = new Hashtable<String,Float>();
        Hashtable<String,Float> rtOfHighestArea = new Hashtable<String,Float>();
        
        //when there are hits where there is no chain information available, the "analyte species" has to be added;
        if (StaticUtils.isThereChainInformationAvailable(molNameWORt,allSets)){
          Hashtable<String,Double> totalAreasOfMolSpecies = new Hashtable<String,Double>();
          Vector<String> modsWoChainAssignment = new Vector<String>();
          for (String mod : adductsSorted.keySet()){
            if (!allMods.containsKey(mod))
              continue;
            if (!StaticUtils.isThereChainInformationAvailable(molNameWORt,allMods.get(mod)))
              modsWoChainAssignment.add(mod);
          }
          for (String mod : adductsSorted.keySet()){
            if (!allMods.containsKey(mod))
              continue;
            Vector<LipidParameterSet> sets = allMods.get(mod);
            double areaOfMod = 0d;
            Hashtable<String,Double> molSpeciesAreasOfMod = new Hashtable<String,Double>();
            Hashtable<String,Float> highestAreas = new Hashtable<String,Float>();
            for (LipidParameterSet set : sets){
              double areaOfOnePeak = set.getArea(maxIsotope);
              areaOfMod += areaOfOnePeak;
              short evidence = StaticUtils.determineEvidenceStateOfHit(set);
              if (!(set instanceof LipidomicsMSnSet) || ((LipidomicsMSnSet)set).getStatus()<=LipidomicsMSnSet.HEAD_GROUP_DETECTED){
                continue;
              }
              LipidomicsMSnSet msn = (LipidomicsMSnSet)set;
              float highestProbeArea = 0f;
              float rtOfHighestProbe = 0f;
              for (CgProbe probe : msn.getIsotopicProbes().get(0)){
                if (probe.Area>highestProbeArea){
                  highestProbeArea = probe.Area;
                  rtOfHighestProbe = probe.Peak;
                }
              }
              for (Object msnNames : msn.getMSnIdentificationNames()){
                String detection = null;
                String detectionWithAmbuiguities = null;
                if (msnNames instanceof Vector){
                  detection = ((Vector<String>)msnNames).get(0);
                  detectionWithAmbuiguities = "";
                  if (((Vector<String>)msnNames).size()>1){
                    for (String name : (Vector<String>)msnNames){
                      if (detectionWithAmbuiguities.length()>0)
                        detectionWithAmbuiguities += ";";
                      detectionWithAmbuiguities += name;
                    }
                  }
                }else if (msnNames instanceof String){
                  detection = (String)msnNames;
                }
                if (detection==null)
                  throw new ExportException("It is not possible that one MSn detection is null!");
                
                String woPosition = StaticUtils.sortFASequenceUnassigned(detection.replaceAll(LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS, LipidomicsConstants.CHAIN_SEPARATOR_NO_POS),
                    LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
                //this is for removal of unassigned positions, which are at the end of the string
                // ths proceeding does not support the export of DG position assignment
                while (woPosition.endsWith("_-"))
                  woPosition = woPosition.substring(0,woPosition.length()-2);
                
                //first I have to check whether this molecular species already exists
                String speciesInList = isSpeciesAlreadyInList(woPosition,molecularSpecies);
                Hashtable<String,Integer> positions = new Hashtable<String,Integer>();
                int count = 0;
                float highestPeakAreaOfAll = 0f;
                float rtOfHighestPeakAreaOfAll = 0f;
                if (speciesInList!=null){
                  woPosition = speciesInList;
                  positions = positionCount.get(woPosition);
                  count = molSpeciesCount.get(woPosition);
                  if (highestPeakArea.containsKey(woPosition)){
                    highestPeakAreaOfAll = highestPeakArea.get(woPosition);
                    rtOfHighestPeakAreaOfAll = rtOfHighestArea.get(woPosition);
                  }
                }else{
                  molecularSpecies.add(woPosition);
                }
                double relativeArea = 0d;
                if (molSpeciesAreasOfMod.containsKey(woPosition))
                  relativeArea = molSpeciesAreasOfMod.get(woPosition);
                double relativePercentage = msn.getRelativeIntensity(detection);
                relativeArea += relativePercentage*areaOfOnePeak;
                molSpeciesAreasOfMod.put(woPosition,relativeArea);
                //for position
                if (adductsSorted.containsKey(mod) && adductsSorted.get(mod)){
                  count++;
                  if (detection.contains("/")){
                    if (detectionWithAmbuiguities!=null){
                      Vector<String> ambiguous = new Vector<String>();
                      if (ambiguousOnes.containsKey(woPosition))
                        ambiguous = ambiguousOnes.get(woPosition);
                      ambiguous.add(detectionWithAmbuiguities);
                      ambiguousOnes.put(woPosition, ambiguous);
                    }else{
                      int posCount = 0;
                      if (positions.containsKey(detection))
                        posCount = positions.get(detection);
                      posCount++;
                      positions.put(detection,posCount);
                    }
                  }
                }
                positionCount.put(woPosition,positions);
                molSpeciesCount.put(woPosition, count);
                if (relativePercentage*highestProbeArea>highestPeakAreaOfAll){
                  highestPeakAreaOfAll = ((float)relativePercentage)*highestProbeArea;
                  rtOfHighestPeakAreaOfAll = rtOfHighestProbe;
                }
                highestPeakArea.put(woPosition,highestPeakAreaOfAll);
                rtOfHighestArea.put(woPosition,rtOfHighestPeakAreaOfAll);
                
                //for detecting the strongest peak and storing its RT separately for every modification
                float highestArea = 0f;
                if (highestAreas.containsKey(woPosition))
                  highestArea = highestAreas.get(woPosition);
                if ((relativePercentage*areaOfOnePeak)>highestArea){
                  highestAreas.put(woPosition, highestArea);
                  if (!rtsOfMods.containsKey(woPosition))
                    rtsOfMods.put(woPosition, new Hashtable<String,Hashtable<String,Double>>());
                  if (!rtsOfMods.get(woPosition).containsKey(mod))
                    rtsOfMods.get(woPosition).put(mod, new Hashtable<String,Double>());
                  rtsOfMods.get(woPosition).get(mod).put(expName, Double.parseDouble(set.getRt()));
                  for (String aMod : modsWoChainAssignment){
                    if (!rtsOfMods.get(woPosition).containsKey(aMod))
                      rtsOfMods.get(woPosition).put(aMod,new Hashtable<String,Double>());
                    if (rtsOfMods.get(woPosition).get(aMod).containsKey(expName) || !allMods.containsKey(aMod) || allMods.get(aMod).size()==0)
                      continue;
                    @SuppressWarnings("rawtypes")
                    Vector rtAndArea = getRtOfHighestZeroIsoArea(allMods.get(aMod));
                    rtsOfMods.get(woPosition).get(aMod).put(expName,((Double)rtAndArea.get(0))/60d);
                  }
                }
                
                Hashtable<String,String> mods = new Hashtable<String,String>();
                if (modsOfSpecies.containsKey(woPosition))
                  mods = modsOfSpecies.get(woPosition);
                mods.put(mod, mod);
                for (String aMod : modsWoChainAssignment)
                  mods.put(aMod, aMod);
                modsOfSpecies.put(woPosition, mods);
                
                storeEvidenceToHash(evidence, evidenceReliabilityOfMods, woPosition, mod, expName);
              }
            }
            //calculate the relative percentage of one mol species and multiply it with the total area of one modification
            double totalMsnModArea = 0d;
            for (Double area : molSpeciesAreasOfMod.values())
              totalMsnModArea += area;
            for (String molSpecies : molSpeciesAreasOfMod.keySet()){
              double percentage = molSpeciesAreasOfMod.get(molSpecies)/totalMsnModArea;
              double totalAreaOfOneMolSpecies = 0d;
              if (totalAreasOfMolSpecies.containsKey(molSpecies))
                totalAreaOfOneMolSpecies = totalAreasOfMolSpecies.get(molSpecies);
              totalAreaOfOneMolSpecies += percentage*areaOfMod;
              totalAreasOfMolSpecies.put(molSpecies, totalAreaOfOneMolSpecies);
            }
          }
          double totalArea = 0d;
          for (double area : totalAreasOfMolSpecies.values())
            totalArea += area;
          //this is the hash with the percentual values
          Hashtable<String,Double> splitsForExp = new Hashtable<String,Double>();
          for (String molSpecies : totalAreasOfMolSpecies.keySet())
            splitsForExp.put(molSpecies, totalAreasOfMolSpecies.get(molSpecies)/totalArea);
          percentalSplits.put(expName, splitsForExp);
        }else{
          if (species.size()==0){
            species.add(molName);
          }
          Hashtable<String,Double> splitsForExp = new Hashtable<String,Double>();
          splitsForExp.put(molName, 1d);
          percentalSplits.put(expName, splitsForExp);

          //for the found modifications and retention time
          Hashtable<String,String> mods = new Hashtable<String,String>();
          if (modsOfSpecies.containsKey(molName))
            mods = modsOfSpecies.get(molName);
          float highestPeakAreaOfAll = 0f;
          float rtOfHighestPeakAreaOfAll = 0f;
          short evidence;
          short currentEvidence;
          for (String aMod : allMods.keySet()){
            mods.put(aMod, aMod);
            evidence = SummaryVO.EVIDENCE_MS1_ONLY;
            for (LipidParameterSet set : allMods.get(aMod)){
              if (mzTabReliabilityOfSumSpecies == 3){
                if (set instanceof LipidomicsMSnSet && ((LipidomicsMSnSet)set).getStatus()>=LipidomicsMSnSet.HEAD_GROUP_DETECTED)
                  mzTabReliabilityOfSumSpecies = 2;
              }
              currentEvidence = StaticUtils.determineEvidenceStateOfHit(set);
              if (currentEvidence>evidence)
            	evidence = currentEvidence;
            }
            storeEvidenceToHash(evidence, evidenceReliabilityOfMods, molName, aMod, expName); 

            @SuppressWarnings("rawtypes")
            Vector rtAndArea = getRtOfHighestZeroIsoArea(allMods.get(aMod));
            double rtOfHighest = ((Double)rtAndArea.get(0))/60d;
            float highestArea = (Float)rtAndArea.get(1);
            if (highestArea>highestPeakAreaOfAll){
              highestPeakAreaOfAll = highestArea;
              rtOfHighestPeakAreaOfAll = (float)(rtOfHighest*60d);
            }
            if (!rtsOfMods.containsKey(molName))
              rtsOfMods.put(molName, new Hashtable<String,Hashtable<String,Double>>());
            if (!rtsOfMods.get(molName).containsKey(aMod))
              rtsOfMods.get(molName).put(aMod, new Hashtable<String,Double>());
            rtsOfMods.get(molName).get(aMod).put(expName, rtOfHighest);
          }
          rtOfHighestArea.put(molName, rtOfHighestPeakAreaOfAll);    
              
          modsOfSpecies.put(molName, mods);

        }
        highestRts.put(expName, rtOfHighestArea);
      }
      
      if (speciesType==LipidomicsConstants.EXPORT_ANALYTE_TYPE_POSITION){
        positionToWithoutHash = replaceSpeciesByPositionAccordingToEvidence(molecularSpecies, molSpeciesCount,
            positionCount, ambiguousOnes);
        molecularSpecies = new Vector<String>();
        for (String molSpecies : positionToWithoutHash.keySet())
          molecularSpecies.add(molSpecies); 
      }
      
      for (String aSpecies : species){
        //the next line is for the export of MSn data only
//        if (rtsOfMods.get(molName)==null)
//          continue;
        Hashtable<String,Double> areas = calculateRelativeAreas(aSpecies,resultsMol,percentalSplits);
        Float rTime = getRtOfHighestPeak(aSpecies,areas,highestRts);
        Vector<String> adducts = getSortedModifications(adductsSorted,modsOfSpecies.get(aSpecies));
        Vector<Integer> featureRefs = getFeatureRefs(aSpecies,adducts,features);
        SummaryVO sumVO = new SummaryVO(currentSummaryId,aSpecies,null,featureRefs,chemFormula,neutralMassTheoretical,rTime,
            adducts,mzTabReliabilityOfSumSpecies,areas,rtsOfMods.get(molName),evidenceReliabilityOfMods.get(molName),expsOfGroup);
        currentSummaryId++;
        speciesSummaries.add(sumVO);
      }
      
      Hashtable<String,SummaryVO> summariesMolSpecies = new Hashtable<String,SummaryVO>();
      List<DoubleStringVO> molSpeciesTotalAreas = new ArrayList<DoubleStringVO>();
      for (String aSpecies : molecularSpecies){
        String speciesInHash = positionToWithoutHash==null ? aSpecies : positionToWithoutHash.get(aSpecies);
        Hashtable<String,Double> areas = calculateRelativeAreas(speciesInHash,resultsMol,percentalSplits);
        Float rTime = getRtOfHighestPeak(speciesInHash,areas,highestRts);
        Vector<String> adducts = getSortedModifications(adductsSorted,modsOfSpecies.get(speciesInHash));
        Vector<Integer> featureRefs = getFeatureRefs(aSpecies,adducts,features);
        String displayString = aSpecies;
        if (exportDoubleBondPositionsForClass) {
          displayString = getDoubleBondPositionDisplayString(aSpecies, allAssignedDoubleBondPositions);
        }
        SummaryVO sumVO = new SummaryVO(molName,displayString,featureRefs,chemFormula,neutralMassTheoretical,rTime,
            adducts,2,areas,rtsOfMods.get(speciesInHash),evidenceReliabilityOfMods.get(speciesInHash),expsOfGroup);
        summariesMolSpecies.put(aSpecies, sumVO);
        double totalArea = 0d;
        for (Double area : areas.values()) totalArea += area;
        molSpeciesTotalAreas.add(new DoubleStringVO(aSpecies,totalArea));
      }
      //sort the areas from the higher to the lower abundant molecular species
      Collections.sort(molSpeciesTotalAreas,new GeneralComparator("at.tugraz.genome.lda.vos.DoubleStringVO", "getValue", "java.lang.Double"));
      for (int i=(molSpeciesTotalAreas.size()-1); i!=-1; i--){
        SummaryVO sumVO = summariesMolSpecies.get(molSpeciesTotalAreas.get(i).getKey());
        sumVO.setId(currentSummaryId);
        currentSummaryId++;
        speciesSummaries.add(sumVO);
      }
    }
    
    //this is for generating the evidence VOs
    Vector<EvidenceVO> evidence = null;
    if (extractEvidence && extractFeatures){
      evidence = new Vector<EvidenceVO>();
      for (String mod : features.keySet()){
        FeatureVO feature = features.get(mod);
        EvidenceBase base = new EvidenceBase(mod,feature.getCharge());
        Vector<String> speciesToCheck = new Vector<String>();
        speciesToCheck.add(molName);
        Hashtable<String,SummaryVO> speciesToSummary = new Hashtable<String,SummaryVO>();
        if (speciesSummaries.size()>0 && speciesSummaries.get(0).getMolecularId()==null)
          speciesToSummary.put(molName, speciesSummaries.get(0));
        for (SummaryVO summary : speciesSummaries){
          if (summary.getMolecularId()!=null){
            speciesToCheck.add(summary.getMolecularId());
            speciesToSummary.put(summary.getMolecularId(), summary);
          }
        }
        //first key: experiment; second key: msLevel; third key: grouping id; value: the scan numbers;
        Hashtable<String,Hashtable<Integer,Hashtable<Integer,Set<Integer>>>> evGroupingIdLookup = new Hashtable<String,Hashtable<Integer,Hashtable<Integer,Set<Integer>>>>();
        for (String species : speciesToCheck){
          String molecularSpecies = null;
          String speciesId;
          String ldaStructure;
          if (speciesToSummary.containsKey(species) && speciesToSummary.get(species).getMolecularId()!=null)
            molecularSpecies = speciesToSummary.get(species).getMolecularId();
          String noPosition = null;
          boolean containsPositionInfo = false;
          if (molecularSpecies!=null){
           noPosition = StaticUtils.sortFASequenceUnassigned(molecularSpecies.replaceAll(LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS, LipidomicsConstants.CHAIN_SEPARATOR_NO_POS),LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
            while (noPosition.endsWith("_-"))
              noPosition = noPosition.substring(0,noPosition.length()-2);
            if (molecularSpecies.indexOf(LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS)!=-1 && !StaticUtils.areAllChainsTheSame(molecularSpecies))
              containsPositionInfo = true;
          }
          LipidParameterSet oneSet = null;
          String chemFormulaInclAdduct = null; 
          neutralMassTheoretical = null;
          for (String expName : expNames){
            if (!evGroupingIdLookup.containsKey(expName))
              evGroupingIdLookup.put(expName, new Hashtable<Integer,Hashtable<Integer,Set<Integer>>>());
            speciesId = new String(molName);
            ldaStructure = null;
            boolean hasHitPositionInfo = false;
            if (!relevantOriginals.containsKey(expName) || !relevantOriginals.get(expName).containsKey(mod) ||
                (speciesToSummary.containsKey(species) && speciesToSummary.get(species).getArea(expName)==0d)){
              continue;
            }
            Vector<LipidParameterSet> sets = relevantOriginals.get(expName).get(mod);
            boolean evaluate = false;
            if (molecularSpecies==null){
              //when the species is the base species, there must not be any molecular identifications present to be reported here
              if (speciesType==LipidomicsConstants.EXPORT_ANALYTE_TYPE_SPECIES || !StaticUtils.isThereChainInformationAvailable(molNameWORt,sets))
                evaluate = true;
            }else{
              if (StaticUtils.isThereChainInformationAvailable(molNameWORt,sets))
                evaluate = true;
            }
            // do not evaluate a species that contains chain information,
            // and do not evaluate a molecular species that does not contain chain information  
            if (!evaluate)
              continue;
            String bestIdentification = null;
            Set<Integer> msLevels = new HashSet<Integer>();
            double weightedMz = 0d;
            double totalArea = 0d;
            double areaTimesMz = 0d;
            boolean foundOneHit = false;
            Hashtable<Integer,Set<Integer>> scanNrs = new Hashtable<Integer,Set<Integer>>();
            
            for (LipidParameterSet set : relevantOriginals.get(expName).get(mod)){
              if (!(set instanceof LipidomicsMSnSet) || ((LipidomicsMSnSet)set).getStatus()<LipidomicsMSnSet.HEAD_GROUP_DETECTED)
                continue;
              evaluate = false;
              oneSet = set;
              LipidomicsMSnSet msn = (LipidomicsMSnSet)set;
              //when the molecular species is found - the set must contain it
              if (molecularSpecies!=null){
                for (Object nameObject : msn.getMSnIdentificationNames()){
                  String detection = null;
                  if (nameObject instanceof Vector){
                    detection = ((Vector<String>)nameObject).get(0);
                  }else if (nameObject instanceof String){
                    detection = (String)nameObject;
                  }
                  if (detection==null)
                    throw new ExportException("It is not possible that one MSn detection is null!");
                  String woPosition = StaticUtils.sortFASequenceUnassigned(detection.replaceAll(LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS, LipidomicsConstants.CHAIN_SEPARATOR_NO_POS),
                      LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
                  //this is for removal of unassigned positions, which are at the end of the string
                  // ths proceeding does not support the export of DG position assignment
                  while (woPosition.endsWith("_-"))
                    woPosition = woPosition.substring(0,woPosition.length()-2);
                  if (StaticUtils.isAPermutedVersion(woPosition,noPosition,LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)){
                    evaluate = true;
                    Set<Integer> msnLevels = msn.getMSLevels(detection,faHydroxyEncoding,lcbHydroxyEncoding);
                    msLevels.addAll(msnLevels);
                    addMsnScanNrsToHash(expName,msnLevels,msn.getMsnRetentionTimes(),scanNrs);
                    //for setting the best identification for found hits of an experiment
                    String identification = detection;
                    if (nameObject instanceof Vector){
                      identification = "";
                      for (String name : (Vector<String>)nameObject){
                        if (identification.length()>0) identification += " | ";
                        identification += name;
                      }
                    }
                    if (bestIdentification==null)
                      bestIdentification = identification;
                    else if (bestIdentification.contains("/")){
                      if (bestIdentification.contains("|") && identification.contains("/") && !identification.contains("|"))
                        bestIdentification = identification;
                    } else if (identification.contains("/"))
                      bestIdentification = identification;

                    
                    if (containsPositionInfo && detection.indexOf("/")!=-1 && hasPositionEvidence(msn.getPositionEvidence(detection,
                        faHydroxyEncoding,lcbHydroxyEncoding)))
                      hasHitPositionInfo = true;
                    break;
                  }
                }         
              }else{
                evaluate = true;
                bestIdentification = molNameWORt;
                if (speciesType==LipidomicsConstants.EXPORT_ANALYTE_TYPE_SPECIES){
                  Set<Integer> msnLevels = msn.getMSLevels(LipidomicsMSnSet.MSLEVEL_ALL_IDENTIFIER,faHydroxyEncoding,lcbHydroxyEncoding);
                  msLevels.addAll(msnLevels);
                  addMsnScanNrsToHash(expName,msnLevels,msn.getMsnRetentionTimes(),scanNrs);
                } else if (!StaticUtils.isThereChainInformationAvailable(molNameWORt,sets)){
                  Set<Integer> msnLevels = msn.getMSLevels(LipidomicsMSnSet.MSLEVEL_HEAD_IDENTIFIER,faHydroxyEncoding,lcbHydroxyEncoding);
                  msLevels.addAll(msnLevels);
                  addMsnScanNrsToHash(expName,msnLevels,msn.getMsnRetentionTimes(),scanNrs);
                }
              }
              if (!evaluate)
                continue;
              foundOneHit = true;
              double areaOfOnePeak = set.getArea(maxIsotope);
              //calculate for each set a weighted mean
              double areaOfSet = 0d;
              double areaTimesMzOfSet = 0d;
              for (CgProbe probe : set.getIsotopicProbes().get(0)){
                areaOfSet += (double)probe.Area;
                areaTimesMzOfSet += (double)probe.Area*(double)probe.Mz;
              }
              double weightedMeanArea = areaTimesMzOfSet/areaOfSet;
              totalArea += areaOfOnePeak;
              areaTimesMz += weightedMeanArea*areaOfOnePeak;
            }
            if (foundOneHit){
              if (chemFormulaInclAdduct==null){
                try {
                  Hashtable<String,Integer> categorized = StaticUtils.categorizeFormula(oneSet.getChemicalFormula());
                  chemFormulaInclAdduct = StaticUtils.getFormulaInHillNotation(categorized,false);
                  neutralMassTheoretical = Settings.getElementParser().calculateTheoreticalMass(StaticUtils.getFormulaInHillNotation(categorized,true), false);
                }catch (ChemicalFormulaException e) {e.printStackTrace();}
              }
              if (molecularSpecies!=null && bestIdentification.indexOf("/")==-1)
                bestIdentification = StaticUtils.sortFASequenceUnassigned(bestIdentification,LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
              
              weightedMz = areaTimesMz/totalArea;
              if (molecularSpecies!=null){
                if (containsPositionInfo && !hasHitPositionInfo)
                  ldaStructure = noPosition;
                else
                  ldaStructure = bestIdentification;
              }

              
              List<Integer> msLevelsSorted = new ArrayList<Integer>(msLevels);              
              for (int msLevel : msLevelsSorted){
                if (!evGroupingIdLookup.get(expName).containsKey(msLevel))
                  evGroupingIdLookup.get(expName).put(msLevel, new Hashtable<Integer,Set<Integer>>());
                int currentGroupingId = checkForSameGroupingEvidence(evGroupingIdLookup.get(expName).get(msLevel),scanNrs.get(msLevel));
                if (currentGroupingId==-1){
                  currentGroupingId = currentEvGroupingId;
                  evGroupingIdLookup.get(expName).get(msLevel).put(currentEvGroupingId, scanNrs.get(msLevel));
                  currentEvGroupingId++;
                }
                evidence.add(new EvidenceVO(expName,base,currentEvidenceId,currentGroupingId,speciesId,ldaStructure,chemFormulaInclAdduct,weightedMz,
                    neutralMassTheoretical,msLevel,scanNrs.get(msLevel)));
                feature.addEvidenceRef(currentEvidenceId);
                currentEvidenceId++;
              }
            }
          }
        }
      }      
    }
       
    SpeciesExportVO exportVO = new SpeciesExportVO(currentSummaryId,speciesSummaries,currentFeatureId,
        features!=null ? new Vector<FeatureVO>(features.values()) : new Vector<FeatureVO>(),currentEvidenceId,currentEvGroupingId,evidence);
    return exportVO;
  }
  
  /**
   * @param aSpecies human readable molecular species of this analyte without double bond position information
   * @param allAssignedDoubleBondPositions assigned double bond positions of this analyte
   * @return Human readable molecular species with double bond position information and the sn position of 'aSpecies'
   */
  public static String getDoubleBondPositionDisplayString(String aSpecies, Vector<DoubleBondPositionVO> allAssignedDoubleBondPositions) {
    String displayString = aSpecies;
    //snPositionLevel adjusts sn position information of the double bond positions display String to the result of replaceSpeciesByPositionAccordingToEvidence
    int snPositionLevel = (aSpecies.contains(LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS)) ?
        LipidomicsConstants.EXPORT_ANALYTE_TYPE_POSITION : LipidomicsConstants.EXPORT_ANALYTE_TYPE_CHAIN;
    String previousAssignedSpecies = null;
    for (DoubleBondPositionVO assignedDoubleBondPosition : allAssignedDoubleBondPositions) {
      String assignedSpecies = assignedDoubleBondPosition.getDoubleBondPositionsHumanReadable(0);
      if (StaticUtils.isChainCombinationEquivalent(assignedSpecies, aSpecies)) {
        if (previousAssignedSpecies != null && !previousAssignedSpecies.equals(assignedSpecies)) {
          System.out.println(String.format("Analytes with differing assigned double bond positions have been grouped together. (%s and %s) "
              + "Consider lowering the setting for 'Show hits with different RT separately' in the 'Selection' tab.", previousAssignedSpecies, assignedSpecies));
          displayString = aSpecies;
          break;                
//          throw new RetentionTimeGroupingException(
//              "<html>Analytes with differing assigned double bond positions have been falsely grouped together.<br/>"
//              + "Consider lowering the setting for 'Show hits with different RT separately' in the 'Selection' tab.</html>");
        }
        if (snPositionLevel == LipidomicsConstants.EXPORT_ANALYTE_TYPE_CHAIN) {
          displayString = assignedSpecies;
        } else if (snPositionLevel == LipidomicsConstants.EXPORT_ANALYTE_TYPE_POSITION && assignedDoubleBondPosition.areChainPositionsFixed()) {
          displayString = assignedDoubleBondPosition.getDoubleBondPositionsHumanReadable(snPositionLevel);
        } else if (snPositionLevel == LipidomicsConstants.EXPORT_ANALYTE_TYPE_POSITION && assignedDoubleBondPosition.getMolecularSpecies().contains(LipidomicsConstants.CHAIN_SEPARATOR_KNOWN_POS)) {
          displayString = assignedSpecies;
        }
        previousAssignedSpecies = assignedSpecies;
      }
    }
    return displayString;
  }

  /**
   * returns a hashtable of adducts sorted by their abundance for each lipid class
   * @param lClasses the available lipid classes;
   * @param originalExcelResults the results from the LDA-resultExcelFiles
   * @return the primary adducts for the lipid classes
   */
  public static Hashtable<String,LinkedHashMap<String,Boolean>> extractAdductsSortedByAbundance(Collection<String> lClasses,
      Hashtable<String,QuantificationResult> originalExcelResults){
    Hashtable<String,LinkedHashMap<String,Boolean>> adductsSorted = new Hashtable<String,LinkedHashMap<String,Boolean>>();
    
    for (String aClass : lClasses){
      LinkedHashMap<String,Boolean> modsSorted = extractAdductsSortedByAbundance(aClass, originalExcelResults);
      if (modsSorted==null)
        continue;
      adductsSorted.put(aClass, modsSorted);
    }
    
    return adductsSorted;
  }
  
  /**
   * returns adducts sorted by their abundance of a specific lipid class
   * @param aClass the lipid class under observation;
   * @param originalExcelResults the results from the LDA-resultExcelFiles
   * @return the primary adducts for the lipid classes
   */
  @SuppressWarnings("unchecked")
  protected static LinkedHashMap<String,Boolean> extractAdductsSortedByAbundance(String aClass,
      Hashtable<String,QuantificationResult> originalExcelResults){
    Hashtable<String,DoubleStringVO> areasPerAdduct = new Hashtable<String,DoubleStringVO>();
    Hashtable<String,Boolean> containsPositionInformation = new Hashtable<String,Boolean>();
    double area;
    for (QuantificationResult result : originalExcelResults.values()){
      if (!result.getIdentifications().containsKey(aClass))
        continue;
      for (LipidParameterSet hit : result.getIdentifications().get(aClass)){
        area = 0d;
        if (areasPerAdduct.containsKey(hit.getModificationName()))
          area = areasPerAdduct.get(hit.getModificationName()).getValue();
        area+=hit.getArea();
        if (area>0){
          areasPerAdduct.put(hit.getModificationName(), new DoubleStringVO(hit.getModificationName(),area));
          boolean posInfo = false;
          if (containsPositionInformation.containsKey(hit.getModificationName()))
            posInfo = containsPositionInformation.get(hit.getModificationName());
          if ((hit instanceof LipidomicsMSnSet) && ((LipidomicsMSnSet)hit).getStatus()>=LipidomicsMSnSet.POSITION_DETECTED &&
              ((LipidomicsMSnSet)hit).getPositionEvidence().size()>0){
           
            boolean isARuleFulfilled = false;             
            for (Hashtable<Integer,Vector<IntensityPositionVO>> posEv : ((LipidomicsMSnSet)hit).getPositionEvidence().values()){
              if (hasPositionEvidence(posEv))
                isARuleFulfilled = true;
              if (isARuleFulfilled)
                break;
            }
            if (isARuleFulfilled)
              posInfo = true;
          } 
          containsPositionInformation.put(hit.getModificationName(), posInfo);
        }
      }
    }
    if (areasPerAdduct.size()==0)
      return null;
    List<DoubleStringVO> unsorted = new ArrayList<DoubleStringVO>(areasPerAdduct.values());
    Collections.sort(unsorted,new GeneralComparator("at.tugraz.genome.lda.vos.DoubleStringVO", "getValue", "java.lang.Double"));
    LinkedHashMap<String,Boolean> modsSorted = new LinkedHashMap<String,Boolean>();
    for (int i=(unsorted.size()-1); i!=-1; i--){
      String adduct = unsorted.get(i).getKey(); 
      boolean positionInformation = false;
      if (containsPositionInformation.containsKey(adduct))
        positionInformation = containsPositionInformation.get(adduct);
      modsSorted.put(adduct,positionInformation);
    }
    return modsSorted;
  }

  
  /**
   * splits the values from the heat maps according to relative molecular species content 
   * @param species the molecular species (without position)
   * @param expResults the results according to the settings in the heat map
   * @param percentualSplits the relative shares how the values from the heat map shall be split
   * @return
   */
  private static Hashtable<String,Double> calculateRelativeAreas(String species,
      Hashtable<String,Vector<Double>> expResults, Hashtable<String,Hashtable<String,Double>> percentualSplits){
    Hashtable<String,Double> relAreas = new Hashtable<String,Double>();
    for (String exp : expResults.keySet()){
      if (expResults.get(exp).get(0)==0d){
        relAreas.put(exp, 0d);
        continue;
      }
      if (percentualSplits!=null){
        if (percentualSplits.get(exp).containsKey(species))
          relAreas.put(exp, expResults.get(exp).get(0)*percentualSplits.get(exp).get(species));
      } else
        relAreas.put(exp, expResults.get(exp).get(0));
    }
    return relAreas;
  }

  /**
   * searches the experiments for a certain species and returns the retention time of the apex of the highest peak
   * @param species the species to look for
   * @param areas the areas to search
   * @param highestRts first key: experiment name; second key: molecular species name (without position); value the retention time of this hit
   * @return
   */
  private static float getRtOfHighestPeak(String species, Hashtable<String,Double> areas, Hashtable<String,Hashtable<String,Float>> highestRts){
    float rt = 0f;
    String strongestExp = getStrongestExp(areas);
    if (strongestExp!=null){
      rt = highestRts.get(strongestExp).get(species);
    }
    return rt;
  }
  
  /**
   * returns the name of the experiment containing the strongest area
   * @param areas key: experiment; value the strongest area in this experiment
   * @return the name of the experiment containing the strongest area
   */
  private static String getStrongestExp(Hashtable<String,Double> areas){
    String strongestExp = null;
    double area = 0d;
    for (String exp : areas.keySet()){
      if (areas.get(exp)>area){
        strongestExp = exp;
        area = areas.get(exp);
      }
    }
    return strongestExp;
  }

  
  /**
   * returns a sorted (by abundance) list of modifications
   * @param available all available modifications sorted by abundance
   * @param found the found hits
   * @return a sorted (by abundance) list of modifications
   */
  private static Vector<String> getSortedModifications(LinkedHashMap<String,Boolean> available, Hashtable<String,String> found){
    Vector<String> sorted = new Vector<String>();
    for (String mod : available.keySet()){
      if (found.containsKey(mod))
        sorted.add(mod);
    }
    return sorted;
  }
  
  
  /**
   * returns a list of identifiers to features belonging to this summary object
   * @param displayName the display name of this analyte
   * @param adducts a list of adducts
   * @param features the available features; key: adduct; value: the feature
   * @return a list of identifiers to features belonging to this summary object
   */
  private static Vector<Integer> getFeatureRefs(String displayName, Vector<String> adducts, LinkedHashMap<String,FeatureVO> features){
    if (features==null)
      return null;
    Vector<Integer> featureRefs = new Vector<Integer>();
    for (String adduct : adducts){
      featureRefs.add(features.get(adduct).getId());
    }
    return featureRefs;
  }

  
  /**
   * checks the position entry of an LipidomicsMSnSet whether there is position information available
   * @param posEv position entry of an LipidomicsMSnSet
   * @return true when there is position evidence present
   */
  private static boolean hasPositionEvidence(Hashtable<Integer,Vector<IntensityPositionVO>> posEv){
    if (posEv==null)
      return false;
    boolean isARuleFulfilled = false;
    for (Vector<IntensityPositionVO> rules : posEv.values()){
      if (rules.size()>0){
        isARuleFulfilled = true;
        break;
      }
    }
    return isARuleFulfilled;
  }

  
  /**
   * adds the found MSn scan numbers to an hash; first key MS-Level; second key: experiment name; value: set of scan numbers
   * @param expName the experiment name
   * @param msnLevels the possible MS-levels
   * @param foundScans first key: MS-level; second key: sorted list of scan numbers; value: retention time
   * @param toAdd the hash table where the scan numbers are added; first key MS-Level; value: set of scan numbers
   */
  private static void addMsnScanNrsToHash(String expName, Set<Integer> msnLevels, Hashtable<Integer,LinkedHashMap<Integer,Float>> foundScans,
      Hashtable<Integer,Set<Integer>> toAdd){
    //for grouping the scanNrs to MS-level and experiment names
    for (Integer msLevel : msnLevels){
      Set<Integer> scans = new HashSet<Integer>();
      if (toAdd.containsKey(msLevel))
        scans = toAdd.get(msLevel);
      for (Integer scanNr : foundScans.get(msLevel).keySet())
        scans.add(scanNr);
      toAdd.put(msLevel, scans);
    }
  }
  
  /**
   * iterates over LipidParameterSets and returns the retention time of the strongest peak (first value in Double), and its area (second value in Float)
   * @param params the collection of LipidParameterSets
   * @return 1) retention time of the strongest peak (Double); 2) its area (Float)
   */
  @SuppressWarnings("rawtypes")
  private static Vector getRtOfHighestZeroIsoArea(Collection<LipidParameterSet> params){
    float highestZeroIsoArea = 0f;
    double rTime = -1d;
    for (LipidParameterSet param : params){
      float areaSplit = 1f;
      if (param.getPercentalSplit()>0)
        areaSplit = param.getPercentalSplit();
      for (CgProbe probe : param.getIsotopicProbes().get(0)){
        if (probe.Area*areaSplit>highestZeroIsoArea){
          highestZeroIsoArea = probe.Area*areaSplit;
          rTime = (double)probe.Peak;
        }
      }
    }
    Vector<Object> timeAndArea = new Vector<Object>();
    timeAndArea.add(rTime);
    timeAndArea.add(highestZeroIsoArea);
    return timeAndArea;
  }
  
  /**
   * adds the evidence add the corresponding position in the hash table (and takes care which evidence is stored, in case of several possibilities) 
   * @param evidence the evidence identifier to store
   * @param hash the hash table; first key: species second key modification; third key: experiment
   * @param molName species name
   * @param mod modification/adduct name
   * @param expName experiment name
   */
  private static void storeEvidenceToHash(short evidence, Hashtable<String,Hashtable<String,Hashtable<String,Short>>> hash,
	  String molName, String mod, String expName){ 
    if (!hash.containsKey(molName))
   	  hash.put(molName, new Hashtable<String,Hashtable<String,Short>>());
    if (!hash.get(molName).containsKey(mod))
      hash.get(molName).put(mod, new Hashtable<String,Short>());
    if (!hash.get(molName).get(mod).containsKey(expName) || evidence>hash.get(molName).get(mod).get(expName))
      hash.get(molName).get(mod).put(expName, evidence);
  }
  
  /**
   * checks whether these spectra belong to an already used mzTab grouping id; if so, the grouping id is returned, otherwise -1
   * @param groupIds the used grouping ids and the refered spectra
   * @param spectra the spectra that are checked
   * @return grouping id if the same spectra were used already, otherwise -1
   */
  private static int checkForSameGroupingEvidence(Hashtable<Integer,Set<Integer>> groupIds, Set<Integer> spectra){
    Set<Integer> groupSpectra;
    for (Integer groupId : groupIds.keySet()){
      groupSpectra = groupIds.get(groupId);
      if (spectra.size()!=groupSpectra.size())
        continue;
      boolean allFound = true;
      for (Integer specNr : spectra){
        if (!groupSpectra.contains(specNr)){
          allFound = false;
          break;
        }
      }
      if (allFound)
        return groupId;
    }
    return -1;
  }
}
