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

package at.tugraz.genome.lda.analysis;

import java.util.Hashtable;
import java.util.Vector;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.analysis.vos.IsoLocationSpaces;
import at.tugraz.genome.lda.analysis.vos.IsoOverlapVO;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.parser.spectrummill.ElementConfigParser;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class IsotopeCorrector
{

  private final static int ISO_MAX_CORRECTION = 5;

  public static Hashtable<String,Vector<LipidParameterSet>> correctIsotopicPattern(ElementConfigParser elementParser, Hashtable<String,Vector<LipidParameterSet>> toCorrect) throws SpectrummillParserException{
    Hashtable<String,LipidParameterSet> correctedParams = new Hashtable<String,LipidParameterSet>();
    Hashtable<String,LipidParameterSet> forCorrection = new Hashtable<String,LipidParameterSet>();
    Hashtable<String,IsoOverlapVO> correctionInfo = new Hashtable<String,IsoOverlapVO>();
    for (String group1 : toCorrect.keySet()){
      Vector<LipidParameterSet> mols1 = toCorrect.get(group1);
      for (LipidParameterSet mol : mols1){
        String analId = uniqueId(group1,mol.getNameIncludingModification());
        IsoOverlapVO overlap = new IsoOverlapVO(analId,mol.getIsotopicProbes().size());
        for (String group2 : toCorrect.keySet()){
          Vector<LipidParameterSet> mols2 = toCorrect.get(group2);
          for (LipidParameterSet compare : mols2){
            if (!group1.equalsIgnoreCase(group2)||!mol.getNameIncludingModification().equalsIgnoreCase(compare.getNameIncludingModification())){
              IsoLocationSpaces overlapLocation = checkIsotopicOverlap(elementParser,mol,compare, group2);
              if (overlapLocation.hasOverlap()){
                overlap.addOverlapLocation(overlapLocation);
              }
            }
          }
        }
        if (overlap.hasOverlap()){
          forCorrection.put(overlap.getId(), mol);
          correctionInfo.put(overlap.getId(), overlap);
//          System.out.println(overlap.getId());
//          overlap.printOverlappingAnalytes();
        }else{
          correctedParams.put(analId,mol);
        }
      }
    }
    int latestSize = forCorrection.size();
    Hashtable<String,String> removedAnalytes = new Hashtable<String,String>();
    while (forCorrection.size()>0){
      Hashtable<String,String> correctedIds = new Hashtable<String,String>();
      for (String id : forCorrection.keySet()){
        LipidParameterSet anal = forCorrection.get(id);
        IsoOverlapVO overlap = correctionInfo.get(id);
        if (overlap.isCorrectionPossible(correctedParams,removedAnalytes)){
          if (overlap.reinitOverlapLocations(anal,correctedParams,ISO_MAX_CORRECTION)){
            LipidParameterSet analCorrected = overlap.makeIsotopicCorrection(anal);
            if (analCorrected.getArea()>0){
              correctedParams.put(id,anal);
            }else{
//              System.out.println("Remoooooooooveeeeeeeeeeeeed: "+id);
              removedAnalytes.put(id, id);
            }
          }else{
            correctedParams.put(id,anal);
          }
          correctedIds.put(id, id);
        }else{
//          System.out.println("The analyte "+id+" cannot be corrected, because it overlaps with an overlap!");
//          overlap.printOverlappingAnalytes();
        }
      }
      for (String id : correctedIds.keySet()){
        forCorrection.remove(id);
        correctionInfo.remove(id);
      }
      if (latestSize==forCorrection.size()){
        //TODO: here I have to do somethin with analytes that cannot be corrected!
        break;
      }
      latestSize=forCorrection.size();
    }
    
    
    //here the corrected areas are transfered back to the data representation required.
    Hashtable<String,Vector<LipidParameterSet>> corrected = new Hashtable<String,Vector<LipidParameterSet>>();
    for (String group : toCorrect.keySet()){
      Vector<LipidParameterSet> mols = toCorrect.get(group);
      Vector<LipidParameterSet> corr = new Vector<LipidParameterSet>();
      for (LipidParameterSet mol : mols){
        String analId = uniqueId(group,mol.getNameIncludingModification());
        if (correctedParams.containsKey(analId)) corr.add(correctedParams.get(analId));
        else if (!removedAnalytes.containsKey(analId)){
          corr.add(mol);
          System.out.println("!!!!!!!!!!!!!! "+group+" "+mol.getNameString()+" was not corrected "+analId);
        }
      }
      corrected.put(group, corr);
    }  
    return corrected;
  }
  
  private static IsoLocationSpaces checkIsotopicOverlap(ElementConfigParser elementParser, LipidParameterSet toCheck, LipidParameterSet overlapping, String overlapGroup) throws SpectrummillParserException{
//    IsoOverlapVO over = new IsoOverlapVO(uniqueId(overlapGroup,overlapping.getNameIncludingModification()),ISO_MAX_CORRECTION);
    IsoLocationSpaces over = new IsoLocationSpaces(uniqueId(overlapGroup,overlapping.getNameIncludingModification()),ISO_MAX_CORRECTION);
    if ((toCheck.Mz[0]+ISO_MAX_CORRECTION*LipidomicsConstants.getNeutronMass()+LipidomicsConstants.getCoarseChromMzTolerance(toCheck.Mz[0]))>overlapping.Mz[0]&&(toCheck.Mz[0]-ISO_MAX_CORRECTION*LipidomicsConstants.getNeutronMass()-LipidomicsConstants.getCoarseChromMzTolerance(toCheck.Mz[0]))<overlapping.Mz[0]){
      boolean negative = false;
      Vector<Vector<Double>> bothDistris = elementParser.calculateChemicalFormulaIntensityDistribution(overlapping.getChemicalFormula(), ISO_MAX_CORRECTION, false);
      Vector<Double> distri = bothDistris.get(0);
      if (bothDistris.size()>1){
        Vector<Double> negDistri = bothDistris.get(1);
        if (StaticUtils.useNegativeDistribution(distri,negDistri)){
          distri = negDistri;
          negative = true;
        }
      }
      // the max range is already checked in the previous if and does not need to be checked here again
      if ((!negative&&toCheck.Mz[0]>(overlapping.Mz[0]+LipidomicsConstants.getCoarseChromMzTolerance(toCheck.Mz[0])))||(negative&&(toCheck.Mz[0]<(overlapping.Mz[0]-LipidomicsConstants.getCoarseChromMzTolerance(toCheck.Mz[0]))))){
        over.setParameterSet(overlapping,distri,negative);
        over.determineOverlappingParts(toCheck);
      }  
    }
    return over;
    
  }

  private static String uniqueId(String overlapGroup, String paramName){
    return overlapGroup+"_-"+paramName;
  }
}
