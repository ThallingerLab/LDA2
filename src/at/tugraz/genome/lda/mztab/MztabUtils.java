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

package at.tugraz.genome.lda.mztab;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Vector;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.analysis.ComparativeAnalysis;
import at.tugraz.genome.lda.exception.ExportException;
import at.tugraz.genome.lda.export.LDAExporter;
import at.tugraz.genome.lda.export.vos.EvidenceVO;
import at.tugraz.genome.lda.export.vos.FeatureVO;
import at.tugraz.genome.lda.export.vos.SpeciesExportVO;
import at.tugraz.genome.lda.export.vos.SummaryVO;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.vos.ResultAreaVO;
import at.tugraz.genome.lda.vos.ResultCompVO;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.utils.Calculator;
import de.isas.mztab1_1.model.MsRun;
import de.isas.mztab1_1.model.OptColumnMapping;
import de.isas.mztab1_1.model.Parameter;
import de.isas.mztab1_1.model.SmallMoleculeEvidence;
import de.isas.mztab1_1.model.SmallMoleculeFeature;
import de.isas.mztab1_1.model.SmallMoleculeSummary;
import de.isas.mztab1_1.model.SpectraRef;

/**
 * Class for generating mzTab specific export information 
 *
 * @author Juergen Hartler
 *
 */
public class MztabUtils extends LDAExporter
{
  
  /**
   * extracts the combined information of several experiments according to selections in the heat maps plus detailed data from LDA result files
   * @param speciesType structural level of data export (lipid species, chain level, position level - for details see LipidomicsConstants.EXPORT_ANALYTE_TYPE)
   * @param currentSummaryId a running unique identifier
   * @param currentFeatureId a running unique identifier
   * @param currentEvidenceId a running unique identifier
   * @param currentEvGroupingId a running identifier for evidence originating from the same spectra
   * @param maxIsotopes the maximum isotope that will be used for the export
   * @param analysisModule an ComparativeAnalysis that calculates all values used for the heat map
   * @param msruns the mzTab MSRun objects; key experiment name; value: MSRun object
   * @param originalExcelResults the LDA results from the Excel files; key: experiment name; value LDA result
   * @param molGroup the analyte class
   * @param molName the name of the species (may contain retention time separated by an "_" at the end)
   * @param resultsMol the values according to the heat map selection
   * @param adductsSorted key: the adducts sorted in consecutive manner starting with the strongest representative; value: contains this adduct position information
   * @param expsOfGroup key: group name; value sorted vector of experiments
   * @return a combined object containing the various mzTab specific information
   * @throws ExportException when there is something wrong
   * @throws SpectrummillParserException when there are elements missing in the elementconfig.xml
   */
  public static SmallMztabMolecule createSmallMztabMolecule(short speciesType, int currentSummaryId, int currentFeatureId, int currentEvidenceId,
      int currentEvGroupingId, int maxIsotopes, ComparativeAnalysis analysisModule, Hashtable<String,MsRun> msruns,
      Hashtable<String,QuantificationResult> originalExcelResults, String molGroup, String molName, Hashtable<String,Vector<Double>> resultsMol,
      LinkedHashMap<String,Boolean> adductsSorted, LinkedHashMap<String,Vector<String>> expsOfGroup)
          throws ExportException, SpectrummillParserException{
    
    Hashtable<String,String> modFormulas = new Hashtable<String,String>();
    Hashtable<String,Integer> modCharges = new Hashtable<String,Integer>();
    Hashtable<String,Double> calcMasses = new Hashtable<String,Double>();
    //this set contains only original data which is relevant for this hit
    Hashtable<String,Hashtable<String,Vector<LipidParameterSet>>> relevantOriginals = new Hashtable<String,Hashtable<String,Vector<LipidParameterSet>>>();
    ResultCompVO comp = analysisModule.getResults().get(molGroup).get(molName).values().iterator().next();
    int isotopes = comp.getAvailableIsotopeNr(maxIsotopes);
    boolean isRtGrouped = analysisModule.isRtGrouped();
    Parameter identificationMethod = new Parameter().name("LipidDataAnalyzer").value(Settings.VERSION);
////    Hashtable<String,Vector<Double>> expMasses = new Hashtable<String,Vector<Double>>();
    for (String expName : analysisModule.getExpNamesInSequence()){
      ResultAreaVO areaVO = analysisModule.getResultAreaVO(molGroup,molName,expName);
      if (areaVO!=null){
        if (areaVO.isAStandard())
          isRtGrouped = false;
        Hashtable<String,Vector<LipidParameterSet>> sets = getRelevantOriginalResults(originalExcelResults.get(expName).getIdentifications().get(molGroup),areaVO);
        relevantOriginals.put(expName, sets);
        for (String mod : adductsSorted.keySet()){
          if (!modFormulas.containsKey(mod)&&areaVO.hasModification(mod)){
            modFormulas.put(mod, areaVO.getModificationFormula(mod));
            modCharges.put(mod, areaVO.getCharge(mod));
            calcMasses.put(mod, areaVO.getTheoreticalMass(mod));
          }
        }
      }
    }
    SpeciesExportVO exportVO =  extractExportableSummaryInformation(speciesType,  true, currentSummaryId, currentFeatureId, true, currentEvidenceId,
        currentEvGroupingId, isRtGrouped, adductsSorted, analysisModule.getExpNamesInSequence(),expsOfGroup,  molName, resultsMol,
        relevantOriginals, isotopes);
//    System.out.println("------------------------------------------");
    //generates the mzTab SmallMoleculeSummary section
    Vector<SmallMoleculeSummary> summaries = new Vector<SmallMoleculeSummary>();
    for (SummaryVO vo : exportVO.getSummaries()){
//      System.out.println(id+": "+vo.getNeutralMass());
      SmallMoleculeSummary summary = new SmallMoleculeSummary();
      summary.setSmlId(String.valueOf(vo.getId()));
      List<String> featureRefs = new ArrayList<String>();
      for (Integer featureRef : vo.getFeatureRefs())
        featureRefs.add(featureRef.toString());
      summary.setSmfIdRefs(featureRefs);
      ArrayList<String> chemFormula = new ArrayList<String>();
      chemFormula.add(vo.getChemFormula());
      summary.setChemicalFormula(chemFormula);
      Vector<Double> neutralMasses = new Vector<Double>();
      neutralMasses.add(vo.getNeutralMass());
      summary.setTheoreticalNeutralMass(neutralMasses);
      summary.setRetentionTime(Calculator.roundDBL((double)vo.getRt(),1));
      List<String> adducts = new ArrayList<String>();
      String mztabAdduct;
      for (String adduct : vo.getModifications()){
        mztabAdduct = LipidomicsConstants.getMzTabAdduct(adduct)!=null ? LipidomicsConstants.getMzTabAdduct(adduct) : adduct;
        adducts.add(mztabAdduct);
      }
      summary.setAdductIons(adducts);
      summary.setReliability(String.valueOf(vo.getMzTabReliability()));
      List<Double> abundanceAssay = new ArrayList<Double>();
      for (String expName : analysisModule.getExpNamesInSequence()){
        abundanceAssay.add(vo.getArea(expName));
//        System.out.println(expName+": "+vo.getArea(expName));
      }
      summary.setAbundanceAssay(abundanceAssay);
      List<Double> abundanceStudyVariable = new ArrayList<Double>();
      List<Double> abundanceCoeffvarStudyVariable = new ArrayList<Double>();
      for (String group : expsOfGroup.keySet()){
        abundanceStudyVariable.add(vo.getMeanArea(group));
        abundanceCoeffvarStudyVariable.add(vo.getCoeffVar(group));
      }
      summary.setAbundanceStudyVariable(abundanceStudyVariable);      
      summary.setAbundanceCoeffvarStudyVariable(abundanceCoeffvarStudyVariable);

      List<OptColumnMapping> optList = new ArrayList<OptColumnMapping>();
      optList.add(new OptColumnMapping().identifier("lipid_species").value(molGroup+" "+vo.getSpeciesId()));
      String bestId = molGroup+" ";
      bestId += vo.getMolecularId()!=null ? vo.getMolecularId() : isRtGrouped ? vo.getSpeciesId().substring(0,vo.getSpeciesId().lastIndexOf("_")) : vo.getSpeciesId();
      optList.add(new OptColumnMapping().identifier("lipid_best_id_level").value(bestId));
      summary.setOpt(optList);
      ////summary.setExpMassToCharge(vo.getNeutralMass());      
      summaries.add(summary);
    }
    
    //generates the mzTab SmallMoleculeSummary section
    Vector<SmallMoleculeFeature> features = new Vector<SmallMoleculeFeature>();
    if (exportVO.getFeatures()!=null){
      for (FeatureVO vo : exportVO.getFeatures()){
        SmallMoleculeFeature feature = new SmallMoleculeFeature();
        feature.setSmfId(String.valueOf(vo.getId()));
        List<String> smeRefs = new ArrayList<String>();
        for (Integer id : vo.getEvidenceRefs())
          smeRefs.add(String.valueOf(id));
        if (smeRefs.size()>0)
          feature.setSmeIdRefs(smeRefs);
        String mztabAdduct = LipidomicsConstants.getMzTabAdduct(vo.getAdduct())!=null ? LipidomicsConstants.getMzTabAdduct(vo.getAdduct()) : vo.getAdduct();
        feature.setAdductIon(mztabAdduct);
        feature.setExpMassToCharge(vo.getExpMz());
        int charge = vo.getCharge();
        if (mztabAdduct.endsWith("-"))
          charge = charge*-1;
        feature.setCharge(charge);
        feature.setRetentionTime(vo.getRtApex());
        feature.setRetentionTimeStart(vo.getRtStart());
        feature.setRetentionTimeEnd(vo.getRtEnd());
        feature.setAbundanceAssay(new ArrayList<Double>(vo.getAreas()));
        
        //TODO: these are dummy values that are set in the meantime to produce an output:
//        feature.setComment(new ArrayList<Comment>());
//        //feature.setIsotopomer(new Parameter());
//        feature.setSmeIdRefAmbiguityCode(1);
        
        features.add(feature);
      }
    }
    
    //generates the mzTab SmallMoleculeEvidence section
    Vector<SmallMoleculeEvidence> evidences = new Vector<SmallMoleculeEvidence>();
    if (exportVO.getEvidence()!=null){
      for (EvidenceVO vo : exportVO.getEvidence()){
        SmallMoleculeEvidence evidence = new SmallMoleculeEvidence();
        evidence.setSmeId(String.valueOf(vo.getId()));
        evidence.setEvidenceUniqueId(String.valueOf(vo.getEvidenceGroupingId()));
        evidence.setChemicalFormula(vo.getChemFormula());
        String mztabAdduct = LipidomicsConstants.getMzTabAdduct(vo.getModification())!=null ? LipidomicsConstants.getMzTabAdduct(vo.getModification()) : vo.getModification();
        evidence.setAdductIon(mztabAdduct);
        evidence.setExpMassToCharge(vo.getExpMz());
        int charge = vo.getCharge();
        if (mztabAdduct.endsWith("-"))
          charge = charge*-1;
        evidence.setCharge(charge);
        evidence.setTheoreticalMassToCharge(vo.getTheorMz());
        evidence.setIdentificationMethod(identificationMethod);
        evidence.setMsLevel(new Parameter().cvLabel("MS").cvAccession("MS").name("ms level").value(String.valueOf(vo.getMsLevel())));
        List<SpectraRef> spectraNrs = new ArrayList<SpectraRef>();
        for (Integer scanNr : vo.getScanNrs()){
          spectraNrs.add(new SpectraRef().msRun(msruns.get(vo.getExpName())).reference("index="+scanNr.toString()));
        }
        evidence.setSpectraRef(spectraNrs);
        
        List<OptColumnMapping> optList = new ArrayList<OptColumnMapping>();
        optList.add(new OptColumnMapping().identifier("lda_identification").value(molGroup+" "+vo.getLdaId()));
        evidence.setOpt(optList);
        
      //TODO: these are dummy values that are set in the meantime to produce an output:
        evidence.setIdConfidenceMeasure(new ArrayList<Double>());
        
        evidences.add(evidence);
      }
    }
    
    SmallMztabMolecule result = new SmallMztabMolecule(exportVO.getCurrentSummaryId(),summaries,exportVO.getCurrentFeatureId(),features,
        exportVO.getCurrentEvidenceId(),exportVO.getCurrentEvGroupingId(),evidences); 
    return result;
  }
  
  
  
  
  
}
