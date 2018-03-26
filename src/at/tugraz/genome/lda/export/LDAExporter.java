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

import uk.ac.ebi.pride.jmztab1_1.utils.errors.MZTabException;
import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.export.vos.EvidenceBase;
import at.tugraz.genome.lda.export.vos.EvidenceVO;
import at.tugraz.genome.lda.export.vos.FeatureVO;
import at.tugraz.genome.lda.export.vos.SpeciesExportVO;
import at.tugraz.genome.lda.export.vos.SummaryVO;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.msn.vos.IntensityPositionVO;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.DoubleStringVO;
import at.tugraz.genome.lda.vos.ResultAreaVO;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.voutils.GeneralComparator;

/**
 * This class is the base for exporting information based on selctions in the heat maps plus extracting detailed
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
   * @return the original information relevant to this one data point (rectangle)
   */
  protected static Hashtable<String,Vector<LipidParameterSet>> getRelevantOriginalResults(Vector<LipidParameterSet> all, ResultAreaVO areaVO){
    Hashtable<String,Vector<LipidParameterSet>> results = new Hashtable<String,Vector<LipidParameterSet>>();
    for (LipidParameterSet set : all){
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
   */
  protected static String isSpeciesAlreadyInList(String species, Vector<String> molecularSpecies){
    for (String other : molecularSpecies){
      if (species.equalsIgnoreCase(other) || StaticUtils.isAPermutedVersion(species, other)){
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
        if (detections>1 && detections>(totalCount/2))
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
  

  
  
  @SuppressWarnings("unchecked")
  /**
   * extracts the combined information of several experiments according to selections in the heat maps plus detailed data from LDA result files
   * @param speciesType structural level of data export (lipid species, chain level, position level - for details see LipidomicsConstants.EXPORT_ANALYTE_TYPE)
   * @param extractFeatures should the features information be included in the exported data
   * @param currentFeatureId a running unique identifier
   * @param extractEvidence should the evidence information information be included in the exported data (only possible when features are included)
   * @param currentSummaryId a running unique identifier
   * @param currentEvidenceId a running unique identifier
   * @param currentEvGroupingId a running identifier for evidence originating from the same spectra
   * @param isRtGrouped is there only one entry for each species, or are the species grouped be retention time
   * @param adductsSorted key: the adducts sorted in consecutive manner starting with the strongest representative; value: contains this adduct position information 
   * @param expNames the sorted experiment names
   * @param expsOfGroup key: group name; value: experiments belonging to this group
   * @param group the analyte class
   * @param molName the name of the molecule (inclusively an RT grouping parameter)
   * @param resultsMol the values according to the heat map selection
   * @param relevantOriginals the relevant results from the LDA Excel file
   * @param maxIsotope the maximum isotope that will be used for the export
   * @return value object containing all exportable information
   * @throws MZTabException when there is something wrong
   * @throws SpectrummillParserException when there are elements missing in the elementconfig.xml
   */
  protected static SpeciesExportVO extractExportableSummaryInformation(short speciesType, boolean extractFeatures, int currentSummaryId,
      int currentFeatureId, boolean extractEvidence, int currentEvidenceId, int currentEvGroupingId, boolean isRtGrouped,
      LinkedHashMap<String,Boolean> adductsSorted, Vector<String> expNames, LinkedHashMap<String,Vector<String>> expsOfGroup, String group,
      String molName, Hashtable<String,Vector<Double>> resultsMol, Hashtable<String,Hashtable<String,Vector<LipidParameterSet>>> relevantOriginals,
      int maxIsotope)
          throws MZTabException, SpectrummillParserException{
    if (extractEvidence && !extractFeatures)
      throw new MZTabException("It is not possible to extract the evidence without extracting the features");
    String chemFormula = "";
    Double neutralMassTheoretical = null;
    //reading information that is equal for all the identifications of one analyte
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
    Vector<SummaryVO> speciesSummaries = new Vector<SummaryVO>();
    String molNameWORt = new String(molName);
////    String rt = null;
    if (isRtGrouped){
////      rt = molName.substring(molName.lastIndexOf("_")+1);
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
            areas.add(0d);
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
          areas.add(areaOfAssay);
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
      int mzTabReliability = 3;
      Hashtable<String,String> foundMods = new Hashtable<String,String>();
      for (Hashtable<String,Vector<LipidParameterSet>> modParams : relevantOriginals.values()){
        for (String mod : modParams.keySet()){
          foundMods.put(mod, mod);
          for (LipidParameterSet set : modParams.get(mod)){
            if (mzTabReliability == 2) break;
            if (set instanceof LipidomicsMSnSet && ((LipidomicsMSnSet)set).getStatus()>=LipidomicsMSnSet.HEAD_GROUP_DETECTED)
              mzTabReliability = 2;
          }
        }
      }
      //detect the Rt of the strongest identification
      float rTime = 0f;
      if (strongestExp!=null){
        float highestArea = 0f;
        for (Vector<LipidParameterSet> params : relevantOriginals.get(strongestExp).values()){
          for (LipidParameterSet param : params){
            for (CgProbe probe : param.getIsotopicProbes().get(0)){
              if (probe.Area>highestArea){
                highestArea = probe.Area;
                rTime = probe.Peak;
              }
            }
          }
        }
      }
      Vector<String> adducts = getSortedModifications(adductsSorted,foundMods);
      Vector<Integer> featureRefs = getFeatureRefs(molName,adducts,features);
      SummaryVO sumVO = new SummaryVO(currentSummaryId,molName,null,featureRefs,chemFormula,neutralMassTheoretical,rTime,
          adducts,mzTabReliability,areas,expsOfGroup);
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
      
            
      //this is for the summary information
      for (String expName : expNames){
        if (!relevantOriginals.containsKey(expName))
          continue;
        Hashtable<String,Vector<LipidParameterSet>> allMods = relevantOriginals.get(expName);
        Vector<LipidParameterSet> allSets = new Vector<LipidParameterSet>();
        for ( Vector<LipidParameterSet> sets : allMods.values())
          allSets.addAll(sets);
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
            for (LipidParameterSet set : sets){
              double areaOfOnePeak = set.getArea(maxIsotope);
              areaOfMod += areaOfOnePeak;
              if (!(set instanceof LipidomicsMSnSet) || ((LipidomicsMSnSet)set).getStatus()<=LipidomicsMSnSet.HEAD_GROUP_DETECTED)
                continue;
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
                  throw new MZTabException("It is not possible that one MSn detection is null!");
                
                String woPosition = StaticUtils.sortFASequenceUnassigned(detection);
                //this is for removal of unassigned positions, which are at the end of the string
                // ths proceeding does not support the export of DG position assignment
                while (woPosition.endsWith("_-"))
                  woPosition = woPosition.substring(0,woPosition.length()-2);
                
                //first I have to check whether this molecular species already exists
                String speciesInList = isSpeciesAlreadyInList(woPosition,molecularSpecies);
                Hashtable<String,Integer> positions = new Hashtable<String,Integer>();
                int count = 0;
                float highestArea = 0f;
                float rtOfHighest = 0f;
                if (speciesInList!=null){
                  woPosition = speciesInList;
                  positions = positionCount.get(woPosition);
                  count = molSpeciesCount.get(woPosition);
                  if (highestPeakArea.containsKey(woPosition)){
                    highestArea = highestPeakArea.get(woPosition);
                    rtOfHighest = rtOfHighestArea.get(woPosition);
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
                if (relativePercentage*highestProbeArea>highestArea){
                  highestArea = ((float)relativePercentage)*highestProbeArea;
                  rtOfHighest = rtOfHighestProbe;
                }
                highestPeakArea.put(woPosition,highestArea);
                rtOfHighestArea.put(woPosition,rtOfHighest);
                Hashtable<String,String> mods = new Hashtable<String,String>();
                if (modsOfSpecies.containsKey(woPosition))
                  mods = modsOfSpecies.get(woPosition);
                mods.put(mod, mod);
                for (String aMod : modsWoChainAssignment)
                  mods.put(aMod, aMod);
                modsOfSpecies.put(woPosition, mods);
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
          //for the retention time
          float highestArea = 0f;
          float rtOfHighest = 0f;
          for (LipidParameterSet set : allSets){
            for (CgProbe probe : set.getIsotopicProbes().get(0)){
              if (probe.Area>highestArea){
                highestArea = probe.Area;
                rtOfHighest = probe.Peak;
              }
            }
            if (set instanceof LipidomicsMSnSet && ((LipidomicsMSnSet)set).getStatus()>=LipidomicsMSnSet.HEAD_GROUP_DETECTED)
              mzTabReliabilityOfSumSpecies = 2;
          }
          rtOfHighestArea.put(molName, rtOfHighest);
          //for the found modifications
          Hashtable<String,String> mods = new Hashtable<String,String>();
          if (modsOfSpecies.containsKey(molName))
            mods = modsOfSpecies.get(molName);
          for (String aMod : allMods.keySet())
            mods.put(aMod, aMod);
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
        Hashtable<String,Double> areas = calculateRelativeAreas(aSpecies,resultsMol,percentalSplits);
        Float rTime = getRtOfHighestPeak(aSpecies,areas,highestRts);
        Vector<String> adducts = getSortedModifications(adductsSorted,modsOfSpecies.get(aSpecies));
        Vector<Integer> featureRefs = getFeatureRefs(aSpecies,adducts,features);
        SummaryVO sumVO = new SummaryVO(currentSummaryId,aSpecies,null,featureRefs,chemFormula,neutralMassTheoretical,rTime,
            adducts,mzTabReliabilityOfSumSpecies,
            areas,expsOfGroup);
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
        SummaryVO sumVO = new SummaryVO(molName,aSpecies,featureRefs,chemFormula,neutralMassTheoretical,rTime,
            adducts,2,areas,expsOfGroup);
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
        EvidenceBase base = new EvidenceBase(currentEvGroupingId,mod,feature.getCharge());
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
        boolean foundMsn = false;
        for (String species : speciesToCheck){
          String molecularSpecies = null;
          String ldaId = new String(group+molName);
          if (speciesToSummary.containsKey(species) && speciesToSummary.get(species).getMolecularId()!=null)
            molecularSpecies = speciesToSummary.get(species).getMolecularId();
          String noPosition = null;
          boolean containsPositionInfo = false;
          if (molecularSpecies!=null){
            noPosition = StaticUtils.sortFASequenceUnassigned(molecularSpecies.replaceAll("/", "_"));
            while (noPosition.endsWith("_-"))
              noPosition = noPosition.substring(0,noPosition.length()-2);
            if (molecularSpecies.indexOf("/")!=-1 && !StaticUtils.areAllChainsTheSame(molecularSpecies))
              containsPositionInfo = true;
          }
          double totalArea = 0d;
          double areaTimesMz = 0d;
          double weightedMz = 0d;
          boolean foundOneHit = false;
          boolean hasHitPositionInfo = false;
          LipidParameterSet oneSet = null;
          Set<Integer> msLevels = new HashSet<Integer>(); 
          //the first key is the ms-Level, the second the experiments
          Hashtable<Integer,Hashtable<String,Set<Integer>>> scanNrs = new Hashtable<Integer,Hashtable<String,Set<Integer>>>();
          Hashtable<String,String> highestLDAStructuralEv = new Hashtable<String,String>();
          for (String expName : expNames){
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
                    throw new MZTabException("It is not possible that one MSn detection is null!");
                  String woPosition = StaticUtils.sortFASequenceUnassigned(detection);
                  //this is for removal of unassigned positions, which are at the end of the string
                  // ths proceeding does not support the export of DG position assignment
                  while (woPosition.endsWith("_-"))
                    woPosition = woPosition.substring(0,woPosition.length()-2);
                  if (StaticUtils.isAPermutedVersion(woPosition,noPosition)){
                    evaluate = true;
                    Set<Integer> msnLevels = msn.getMSLevels(detection);
                    msLevels.addAll(msnLevels);
                    addMsnScanNrsToHash(expName,msnLevels,msn.getMsnRetentionTimes(),scanNrs);
                    //for setting the best identification for found hits of an experiment
                    String identification = detection;
                    if (nameObject instanceof Vector){
                      identification = "";
                      for (String name : (Vector<String>)nameObject){
                        if (identification.length()>0) identification += "|";
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

                    
                    if (containsPositionInfo && detection.indexOf("/")!=-1 && hasPositionEvidence(msn.getPositionEvidence(detection)))
                      hasHitPositionInfo = true;
                    break;
                  }
                }         
              }else{
                evaluate = true;
                bestIdentification = molNameWORt;
                if (speciesType==LipidomicsConstants.EXPORT_ANALYTE_TYPE_SPECIES){
                  Set<Integer> msnLevels = msn.getMSLevels(LipidomicsMSnSet.MSLEVEL_ALL_IDENTIFIER);
                  msLevels.addAll(msnLevels);
                  addMsnScanNrsToHash(expName,msnLevels,msn.getMsnRetentionTimes(),scanNrs);
                } else if (!StaticUtils.isThereChainInformationAvailable(molNameWORt,sets)){
                  Set<Integer> msnLevels = msn.getMSLevels(LipidomicsMSnSet.MSLEVEL_HEAD_IDENTIFIER);
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
            if (bestIdentification!=null)
              highestLDAStructuralEv.put(expName, bestIdentification);
          }
          if (foundOneHit){
            if (molecularSpecies!=null){
              if (containsPositionInfo && !hasHitPositionInfo)
                ldaId += "|"+noPosition;
              else
                ldaId += "|"+molecularSpecies;
            }
            weightedMz = areaTimesMz/totalArea;          
            String chemFormulaInclAdduct = null; 
            try {
              Hashtable<String,Integer> categorized = StaticUtils.categorizeFormula(oneSet.getChemicalFormula());
              chemFormulaInclAdduct = StaticUtils.getFormulaInHillNotation(categorized,false);
              neutralMassTheoretical = Settings.getElementParser().calculateTheoreticalMass(StaticUtils.getFormulaInHillNotation(categorized,true), false);
            }catch (ChemicalFormulaException e) {e.printStackTrace();}
            List<Integer> msLevelsSorted = new ArrayList<Integer>(msLevels);
            Collections.sort(msLevelsSorted);
            for (int msLevel : msLevelsSorted){
              evidence.add(new EvidenceVO(base,currentEvidenceId,ldaId,highestLDAStructuralEv,chemFormulaInclAdduct,weightedMz,neutralMassTheoretical,msLevel,scanNrs.get(msLevel)));
              feature.addEvidenceRef(currentEvidenceId);
              currentEvidenceId++;
              foundMsn = true;
            }
          }

        }
        if (foundMsn)
          currentEvGroupingId++;
      }      
    }
       
    SpeciesExportVO exportVO = new SpeciesExportVO(currentSummaryId,speciesSummaries,currentFeatureId,
        new Vector<FeatureVO>(features.values()),currentEvidenceId,currentEvGroupingId,evidence);
    return exportVO;
  }

  /**
   * returns a hashtable of adducts sorted by their abundance for each lipid class
   * @param lClasses the available lipid classes;
   * @param originalExcelResults the results from the LDA-resultExcelFiles
   * @return the primary adducts for the lipid classes
   */
  @SuppressWarnings("unchecked")
  public static Hashtable<String,LinkedHashMap<String,Boolean>> extractAdductsSortedByAbundance(Collection<String> lClasses,
      Hashtable<String,QuantificationResult> originalExcelResults){
    Hashtable<String,LinkedHashMap<String,Boolean>> adductsSorted = new Hashtable<String,LinkedHashMap<String,Boolean>>();
    double area;
    for (String aClass : lClasses){
      Hashtable<String,DoubleStringVO> areasPerAdduct = new Hashtable<String,DoubleStringVO>();
      Hashtable<String,Boolean> containsPositionInformation = new Hashtable<String,Boolean>();
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
        continue;
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
      adductsSorted.put(aClass, modsSorted);
    }
    
    return adductsSorted;
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
   * @param toAdd the hash table where the scan numbers are added; first key MS-Level; second key: experiment name; value: set of scan numbers
   */
  private static void addMsnScanNrsToHash(String expName, Set<Integer> msnLevels, Hashtable<Integer,LinkedHashMap<Integer,Float>> foundScans,
      Hashtable<Integer,Hashtable<String,Set<Integer>>> toAdd){
    //for grouping the scanNrs to MS-level and experiment names
    for (Integer msLevel : msnLevels){
      Hashtable<String,Set<Integer>> scanNrsEachScan = new Hashtable<String,Set<Integer>>();
      if (toAdd.containsKey(msLevel))
        scanNrsEachScan = toAdd.get(msLevel);
      Set<Integer> scans = new HashSet<Integer>();
      if (scanNrsEachScan.containsKey(expName))
        scans = scanNrsEachScan.get(expName);
      for (Integer scanNr : foundScans.get(msLevel).keySet())
        scans.add(scanNr);
      scanNrsEachScan.put(expName, scans);
      toAdd.put(msLevel, scanNrsEachScan);                      
    }
  }
  
}