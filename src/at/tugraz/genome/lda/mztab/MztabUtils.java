/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.exception.RetentionTimeGroupingException;
import at.tugraz.genome.lda.export.LDAExporter;
import at.tugraz.genome.lda.export.vos.EvidenceVO;
import at.tugraz.genome.lda.export.vos.FeatureVO;
import at.tugraz.genome.lda.export.vos.SpeciesExportVO;
import at.tugraz.genome.lda.export.vos.SummaryVO;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.vos.ResultAreaVO;
import at.tugraz.genome.lda.vos.ResultCompVO;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import de.isas.mztab2.model.MsRun;
import de.isas.mztab2.model.OptColumnMapping;
import de.isas.mztab2.model.Parameter;
import de.isas.mztab2.model.SmallMoleculeEvidence;
import de.isas.mztab2.model.SmallMoleculeFeature;
import de.isas.mztab2.model.SmallMoleculeSummary;
import de.isas.mztab2.model.SpectraRef;

/**
 * Class for generating mzTab specific export information 
 *
 * @author Juergen Hartler
 *
 */
public class MztabUtils extends LDAExporter
{
  
  /** the value for the best_id_confidence_measure is always the same*/
  private final static Parameter BEST_CONFID_MEASURE = new Parameter().cvLabel("MS").cvAccession("MS:1002890").name("fragmentation score");
  
  /**
   * extracts the combined information of several experiments according to selections in the heat maps plus detailed data from LDA result files
   * @param speciesType structural level of data export (lipid species, chain level, position level - for details see LipidomicsConstants.EXPORT_ANALYTE_TYPE)
   * @param exportDoubleBondPositionsForMolecule true when double bond positions shall be exported for this lipid class
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
   * @param faHydroxyEncoding the OH encodings of the FA moiety
   * @param lcbHydroxyEncoding the OH encodings of the LCB moiety
   * @return a combined object containing the various mzTab specific information
   * @throws ExportException when there is something wrong
   * @throws SpectrummillParserException when there are elements missing in the elementconfig.xml
   * @throws LipidCombinameEncodingException thrown when a lipid combi id (containing type and OH number) cannot be decoded
   */
  public static SmallMztabMolecule createSmallMztabMolecule(short speciesType, boolean exportDoubleBondPositionsForClass, int currentSummaryId, int currentFeatureId, int currentEvidenceId,
      int currentEvGroupingId, int maxIsotopes, ComparativeAnalysis analysisModule, Hashtable<String,MsRun> msruns,
      Hashtable<String,QuantificationResult> originalExcelResults, String molGroup, String molName, Hashtable<String,Vector<Double>> resultsMol,
      LinkedHashMap<String,Boolean> adductsSorted, LinkedHashMap<String,Vector<String>> expsOfGroup, HydroxyEncoding faHydroxyEncoding,
      HydroxyEncoding lcbHydroxyEncoding)
          throws ExportException, SpectrummillParserException, LipidCombinameEncodingException, RetentionTimeGroupingException{
    
    Hashtable<String,String> modFormulas = new Hashtable<String,String>();
    Hashtable<String,Integer> modCharges = new Hashtable<String,Integer>();
    Hashtable<String,Double> calcMasses = new Hashtable<String,Double>();
    //this set contains only original data which is relevant for this hit
    Hashtable<String,Hashtable<String,Vector<LipidParameterSet>>> relevantOriginals = new Hashtable<String,Hashtable<String,Vector<LipidParameterSet>>>();
    ResultCompVO comp = analysisModule.getResults().get(molGroup).get(molName).values().iterator().next();
    int isotopes = comp.getAvailableIsotopeNr(maxIsotopes);
    boolean isRtGrouped = analysisModule.isRtGrouped();
    Parameter identificationMethod = new Parameter().name("LipidDataAnalyzer").value(Settings.VERSION);
    for (String expName : analysisModule.getExpNamesInSequence()){
      ResultAreaVO areaVO = analysisModule.getResultAreaVO(molGroup,molName,expName);
      if (areaVO!=null){
        if (areaVO.isAStandard())
          isRtGrouped = false;
        Hashtable<String,Vector<LipidParameterSet>> sets = getRelevantOriginalResults(originalExcelResults.get(expName).getIdentifications().get(molGroup),areaVO,false);
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
    SpeciesExportVO exportVO = LDAExporter.extractExportableSummaryInformation(speciesType, exportDoubleBondPositionsForClass, true, currentSummaryId, currentFeatureId, true, currentEvidenceId,
        currentEvGroupingId, isRtGrouped, adductsSorted, analysisModule.getExpNamesInSequence(),expsOfGroup,  molName, resultsMol,
        relevantOriginals, isotopes, faHydroxyEncoding, lcbHydroxyEncoding);
//    System.out.println("------------------------------------------");
    //generates the mzTab SmallMoleculeSummary section
    Double smlArea;
    Vector<SmallMoleculeSummary> summaries = new Vector<SmallMoleculeSummary>();
    Vector<String> dbIdentifiers = null;
    Vector<String> chemicalNames = null;
    for (SummaryVO vo : exportVO.getSummaries()){
//      System.out.println(id+": "+vo.getNeutralMass());
      SmallMoleculeSummary summary = new SmallMoleculeSummary();
      summary.setSmlId(vo.getId());
      summary.setSmfIdRefs(vo.getFeatureRefs());
      dbIdentifiers = new Vector<String>();
      dbIdentifiers.add("lda2:"+molGroup+" "+(isRtGrouped ? vo.getSpeciesId().substring(0,vo.getSpeciesId().lastIndexOf("_")) : vo.getSpeciesId()));
      summary.setDatabaseIdentifier(dbIdentifiers);
      ArrayList<String> chemFormula = new ArrayList<String>();
      chemFormula.add(vo.getChemFormula());
      summary.setChemicalFormula(chemFormula);
      
      chemicalNames = new Vector<String>();
      String bestId = molGroup+" ";
      bestId += vo.getMolecularId()!=null ? vo.getMolecularId() : isRtGrouped ? vo.getSpeciesId().substring(0,vo.getSpeciesId().lastIndexOf("_")) : vo.getSpeciesId();
      chemicalNames.add(bestId);
      summary.setChemicalName(chemicalNames);
      
      Vector<Double> neutralMasses = new Vector<Double>();
      neutralMasses.add(vo.getNeutralMass());
      summary.setTheoreticalNeutralMass(neutralMasses);
      List<String> adducts = new ArrayList<String>();
      String mztabAdduct;
      for (String adduct : vo.getModifications()){
        mztabAdduct = LipidomicsConstants.getMzTabAdduct(adduct)!=null ? LipidomicsConstants.getMzTabAdduct(adduct) : adduct;
        adducts.add(mztabAdduct);
      }
      summary.setAdductIons(adducts);
      summary.setReliability(String.valueOf(vo.getMzTabReliability()));
      summary.setBestIdConfidenceMeasure(BEST_CONFID_MEASURE);
      List<Double> abundanceAssay = new ArrayList<Double>();
      for (String expName : analysisModule.getExpNamesInSequence()){
        smlArea = vo.getArea(expName);
        abundanceAssay.add(smlArea>0d ? vo.getArea(expName) : null);
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
      summary.setAbundanceVariationStudyVariable(abundanceCoeffvarStudyVariable);

      List<OptColumnMapping> optList = new ArrayList<OptColumnMapping>();
      optList.add(new OptColumnMapping().identifier("global_lipid_species").value(molGroup+" "+(isRtGrouped ? vo.getSpeciesId().substring(0,vo.getSpeciesId().lastIndexOf("_")) : vo.getSpeciesId())));
      if (isRtGrouped)
        optList.add(new OptColumnMapping().identifier("global_lipid_lda_species").value(molGroup+" "+vo.getSpeciesId()));
      summary.setOpt(optList);
      summaries.add(summary);
    }
    
    //generates the mzTab SmallMoleculeSummary section
    Vector<SmallMoleculeFeature> features = new Vector<SmallMoleculeFeature>();
    Hashtable<String,Short> polarities = new Hashtable<String,Short>();
    for (String expName : analysisModule.getExpNamesInSequence())
      polarities.put(expName,  SmallMztabMolecule.POLARITY_UNKNOWN);
    if (exportVO.getFeatures()!=null){
      for (FeatureVO vo : exportVO.getFeatures()){
        SmallMoleculeFeature feature = new SmallMoleculeFeature();
        feature.setSmfId(vo.getId());
        if (vo.getEvidenceRefs().size()>0)
          feature.setSmeIdRefs(vo.getEvidenceRefs());
        if (vo.getEvidenceRefs().size()>1)
          feature.setSmeIdRefAmbiguityCode(2);
        String mztabAdduct = LipidomicsConstants.getMzTabAdduct(vo.getAdduct())!=null ? LipidomicsConstants.getMzTabAdduct(vo.getAdduct()) : vo.getAdduct();
        feature.setAdductIon(mztabAdduct);
        feature.setExpMassToCharge(vo.getExpMz());
        int charge = vo.getCharge();
        if (mztabAdduct.endsWith("-"))
          charge = charge*-1;
        feature.setCharge(charge);
        feature.setRetentionTimeInSeconds(vo.getRtApex());
        feature.setRetentionTimeInSecondsStart(vo.getRtStart());
        feature.setRetentionTimeInSecondsEnd(vo.getRtEnd());
        feature.setAbundanceAssay(new ArrayList<Double>(vo.getAreas()));
        for (int i=0; i!=analysisModule.getExpNamesInSequence().size(); i++){
          short polarity = SmallMztabMolecule.POLARITY_UNKNOWN;
          if (vo.getAreas().get(i)==null) continue;
          String expName = analysisModule.getExpNamesInSequence().get(i);
          if (mztabAdduct.endsWith("+"))
            polarity = polarities.get(expName)==SmallMztabMolecule.POLARITY_NEGATIVE ? SmallMztabMolecule.POLARITY_BOTH : SmallMztabMolecule.POLARITY_POSITIVE;
          else if (mztabAdduct.endsWith("-"))
            polarity = polarities.get(expName)==SmallMztabMolecule.POLARITY_POSITIVE ? SmallMztabMolecule.POLARITY_BOTH : SmallMztabMolecule.POLARITY_NEGATIVE;
          polarities.put(expName, polarity);
        }
        //TODO: these are dummy values that are set in the meantime to produce an output:
//        feature.setComment(new ArrayList<Comment>());
//        //feature.setIsotopomer(new Parameter());
//        feature.setSmeIdRefAmbiguityCode(1);
        
        features.add(feature);
      }
    }
    
    //generates the mzTab SmallMoleculeEvidence section
    Vector<SmallMoleculeEvidence> evidences = new Vector<SmallMoleculeEvidence>();
    String dbId = null;
    if (exportVO.getEvidence()!=null){
      for (EvidenceVO vo : exportVO.getEvidence()){
        SmallMoleculeEvidence evidence = new SmallMoleculeEvidence();
        evidence.setSmeId(vo.getId());
        evidence.setEvidenceInputId(String.valueOf(vo.getEvidenceGroupingId()));
        String speciesName = molGroup+" "+(isRtGrouped ? vo.getSpeciesId().substring(0,vo.getSpeciesId().lastIndexOf("_")) : vo.getSpeciesId());
        dbId = "lda2:"+speciesName;
        evidence.setDatabaseIdentifier(dbId);
        evidence.setChemicalFormula(vo.getChemFormula());
        evidence.setChemicalName(vo.getLdaStructure()==null ? speciesName : molGroup+" "+vo.getLdaStructure().replaceAll(" \\| ", " | "+molGroup+" "));
        String mztabAdduct = LipidomicsConstants.getMzTabAdduct(vo.getModification())!=null ? LipidomicsConstants.getMzTabAdduct(vo.getModification()) : vo.getModification();
        evidence.setAdductIon(mztabAdduct);
        evidence.setExpMassToCharge(vo.getExpMz());
        int charge = vo.getCharge();
        if (mztabAdduct.endsWith("-"))
          charge = charge*-1;
        evidence.setCharge(charge);
        evidence.setTheoreticalMassToCharge(vo.getTheorMz());
        evidence.setIdentificationMethod(identificationMethod);
        evidence.setMsLevel(new Parameter().cvLabel("MS").cvAccession("MS:1000511").name("ms level").value(String.valueOf(vo.getMsLevel())));
        List<SpectraRef> spectraNrs = new ArrayList<SpectraRef>();
        for (Integer scanNr : vo.getScanNrs()){
          spectraNrs.add(new SpectraRef().msRun(msruns.get(vo.getExpName())).reference("scan="+scanNr.toString()));
        }
        evidence.setSpectraRef(spectraNrs);
        List<Double> confidenceMeasures = new ArrayList<Double>();
        confidenceMeasures.add(null);
        evidence.setIdConfidenceMeasure(confidenceMeasures);
        evidence.setRank(1);
        
        List<OptColumnMapping> optList = new ArrayList<OptColumnMapping>();
        optList.add(new OptColumnMapping().identifier("global_lipid_species").value(molGroup+" "+(isRtGrouped ? vo.getSpeciesId().substring(0,vo.getSpeciesId().lastIndexOf("_")) : vo.getSpeciesId())));
        if (isRtGrouped)
          optList.add(new OptColumnMapping().identifier("global_lipid_lda_species").value(molGroup+" "+vo.getSpeciesId()));
//        if (speciesType >= LipidomicsConstants.EXPORT_ANALYTE_TYPE_CHAIN)
//          optList.add(new OptColumnMapping().identifier("global_lipid_molecular_species").value(vo.getLdaStructure()==null ? null : molGroup+" "+vo.getLdaStructure().replaceAll(" \\| ", " | "+molGroup+" ")));
        evidence.setOpt(optList);
        evidences.add(evidence);
      }
    }
    
    SmallMztabMolecule result = new SmallMztabMolecule(exportVO.getCurrentSummaryId(),polarities,summaries,exportVO.getCurrentFeatureId(),
        features,exportVO.getCurrentEvidenceId(),exportVO.getCurrentEvGroupingId(),evidences); 
    return result;
  }
  
  
  
  
  
}
