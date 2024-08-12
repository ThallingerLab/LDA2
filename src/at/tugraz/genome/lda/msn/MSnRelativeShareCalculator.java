/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2021 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Vector;

import at.tugraz.genome.lda.exception.LMException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.msn.vos.MSnDebugVO;
import at.tugraz.genome.lda.utils.LMLinearFaCombination;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.FloatStringVO;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.util.FloatMatrix;
import at.tugraz.genome.voutils.GeneralComparator;

/**
 * For splitting MS1 intensities according to MSn fragments
 * @author Juergen Hartler
 *
 */
public class MSnRelativeShareCalculator
{
  
  private Hashtable<String,Vector<FattyAcidVO>> combiFAs_;
  private Hashtable<String,String> combiLookup_;
  private LinkedHashMap<FattyAcidVO,Double> sortedFAs_;
  private Vector<String> sortedCombis_;
  /** the cutoff value relative to the highest chain combination*/
  private double relativeChainCutoff_;
  /** the final result; i.e. the relative intensity split according to the prediction - a kind of percentual value*/
  private Hashtable<String,Double> relativeIntensities_;
  /** true when discards by the relative chain cutoff shall be stored in a debug VO*/
  private boolean debug_;
  /** value object for debug purposes*/
  private MSnDebugVO debugVO_;
  
  /**
   * constructor for predicting relative shares of chain combinations sharing the same fragments
   * @param combiFAs hash table containing the individual chain names; key: name of the chain combination; value: vector containing the individual chains
   * @param chainFragments the detected chain fragments; first key: chain combination name; second key: fragment name; value: the detected fragments and areas
   * @param relativeChainCutoff cutoff value relative to the highest chain combination
   * @param debug true when discards by the relative chain cutoff shall be stored in a debug VO
   * @param debugVO value object for debug purposes
   * @throws LipidCombinameEncodingException thrown when there is something wrong with the lipid name
   */
  public MSnRelativeShareCalculator(Hashtable<String,Vector<FattyAcidVO>> combiFAs, Hashtable<String,Hashtable<String,CgProbe>> chainFragments,
      double relativeChainCutoff, boolean debug, MSnDebugVO debugVO) throws LipidCombinameEncodingException{
    this.combiFAs_ = combiFAs;
    Vector<FattyAcidVO> fas = new Vector<FattyAcidVO>();
    for (String chainId : chainFragments.keySet()) fas.add(StaticUtils.decodeLipidNameForCreatingCombis(chainId));
    sortedFAs_ = new LinkedHashMap<FattyAcidVO,Double>();
    for (FattyAcidVO fa : StaticUtils.sortChainVOs(fas)){
      double area = 0d;
      for (CgProbe probe : chainFragments.get(fa.getChainId()).values()){
        area += (double)probe.Area;
      }
    	  sortedFAs_.put(fa, area);
//      System.out.println("original:\t"+fa+"\t"+area);
    }
    combiLookup_ = new Hashtable<String,String>();
    for (String combi : combiFAs.keySet()){
      combiLookup_.put(StaticUtils.encodeLipidCombi(StaticUtils.sortChainVOs(StaticUtils.decodeLipidNamesFromChainCombi(combi))),combi);
    }
    ////sortedCombis_ = sortCombis(sortedFAs_.keySet(),combiLookup_.keySet(),"");
    ////sortedCombis_ = StaticUtils.sortChainCombinations(new HashSet<String>(combiLookup_.values()));
    sortedCombis_ = StaticUtils.sortChainCombinations(new HashSet<String>(combiLookup_.keySet()));
    this.relativeChainCutoff_ = relativeChainCutoff;
    this.debug_ = debug;
    this.debugVO_ = debugVO;
//    int combiCount = 0;
//    for (String combi : sortedCombis_){
//      System.out.println(combiCount+": "+combi);
//      combiCount++;
//    }
  }
  
  
  /**
   * splits the intensities of the species after LM optimization was performed
   */
  @SuppressWarnings("unchecked")
  public void splitIntensities(){
    int negativeParams = 1;
    List<FloatStringVO> negatives;
    Hashtable<String,Integer> lookup;
    Vector<Integer> toRemove;
    Vector<FattyAcidVO> fasToRemove;
    //int itNumber = 1;
    FloatMatrix params = null;
    LMLinearFaCombination optimizer = null;
    while (negativeParams>0){
      //System.out.println("Iteration: "+itNumber);
      negativeParams = 0;
      optimizer = parameterPrediction();
      params = optimizer.getResultParams();
      negatives = new ArrayList<FloatStringVO>();
      lookup = new Hashtable<String,Integer>();
      float highest = 0f;
      for (int i=0;i!=params.A.length; i++) {
        if (params.A[i][0]>highest)
          highest = params.A[i][0];
      }
      for (int i=0;i!=params.A.length; i++){
        //System.out.println(/*i+"."+*/sortedCombis_.get(i)+"\t\t"+params.A[i][0]);
        if (params.A[i][0]>=relativeChainCutoff_*highest)
          continue;
        if (debug_) debugVO_.addViolatedCombinations(combiLookup_.get(sortedCombis_.get(i)), MSnDebugVO.COMBINATION_LOWER_CHAIN_CUTOFF);
        negativeParams++;
        negatives.add(new FloatStringVO(sortedCombis_.get(i),params.A[i][0]));
        lookup.put(sortedCombis_.get(i), i);
      }
      Collections.sort(negatives,new GeneralComparator("at.tugraz.genome.lda.vos.FloatStringVO", "getValue", "java.lang.Float"));
      toRemove = new Vector<Integer>();
      int maxRemove = negatives.size()/2+negatives.size()%2;
      for (int i=0; i!=maxRemove; i++){
        toRemove.add(lookup.get(negatives.get(i).getKey()));
      }
      Collections.sort(toRemove);
      for (int i=(toRemove.size()-1); i!=-1; i--){
        sortedCombis_.remove(toRemove.get(i).intValue());
      }
      fasToRemove = new Vector<FattyAcidVO>();
      for (FattyAcidVO fa : sortedFAs_.keySet()){
        boolean faFound = false;
        for (String combi : sortedCombis_){
          for (FattyAcidVO otherFA : combiFAs_.get(combiLookup_.get(combi))){
            if (fa.getChainId().equalsIgnoreCase(otherFA.getChainId())){
              faFound = true;
              break;
            }            
          }
          if (faFound)
            break;
        }
        if (!faFound){
          //System.out.println("Removing FA: "+fa);
          fasToRemove.add(fa);
        }
      }
      for (FattyAcidVO fa : fasToRemove){
        sortedFAs_.remove(fa);
      }

      
//      for (FloatStringVO neg : negatives){
//        System.out.println(neg.getKey()+": "+neg.getValue());
//      }
      
      //itNumber++;
      //System.out.println("----------------------------------------------------------");
      //this is for stopping the loop
      //negativeParams = 0;
    }
    float total = 0f;
    for (int i=0;i!=params.A.length; i++) total += params.A[i][0];
    relativeIntensities_ = new Hashtable<String,Double>();
    for (int i=0;i!=params.A.length; i++){
      relativeIntensities_.put(combiLookup_.get(sortedCombis_.get(i)), ((double)params.A[i][0])/total);
    }
//    System.out.println("----------------------------------------------------------");
    
    //for testing the output
//    float[][] values = new float[sortedFAs_.size()][sortedCombis_.size()];
//    int rowCount = 0;
//    for (String fa : sortedFAs_.keySet()){
//      int columnCount = 0;
//      for (String combi : sortedCombis_){
//        int faCount = 0;
//        for (String otherFA : combiFAs_.get(combiLookup_.get(combi))){
//          if (fa.equalsIgnoreCase(otherFA))
//            faCount++;
//        }
//        values[rowCount][columnCount] = faCount;
//        columnCount++;
//      }
//      rowCount++;
//    }
//    int faCount = 0;
//    for (String fa : sortedFAs_.keySet()){
//      try {
//        System.out.println(fa+"\t"+optimizer.calculateFitValue(values[faCount]));
//      }
//      catch (LMException e) {
//        // TODO Auto-generated catch block
//        e.printStackTrace();
//      }
//      faCount++;
//    }
  }
  
  /**
   * applies the Levenberg-Marquardt algorithm to predict the contributions of each individiual fatty acid
   * @return the fitted Levenberg-Marquardt model containing the final parameters - getResultParams() returns the fitted parameters
   */
  private LMLinearFaCombination parameterPrediction(){
    //create the input matrices
    float[][] values = new float[sortedFAs_.size()][sortedCombis_.size()];
    float[] intensities = new float[sortedFAs_.size()];
    int rowCount = 0;
    for (FattyAcidVO fa : sortedFAs_.keySet()){
//      System.out.println("!!! "+fa);
      int columnCount = 0;
      for (String combi : sortedCombis_){
        int faCount = 0;
        for (FattyAcidVO otherFA : combiFAs_.get(combiLookup_.get(combi))){
          if (fa.getChainId().equalsIgnoreCase(otherFA.getChainId()))
            faCount++;
        }
        values[rowCount][columnCount] = faCount;
        columnCount++;
      }
      intensities[rowCount] = sortedFAs_.get(fa).floatValue();
      rowCount++;
    }

    LMLinearFaCombination optimizer = new LMLinearFaCombination(values,intensities);
    try {
      optimizer.fit();
      
      return optimizer;
    }
    catch (LMException e) {
      e.printStackTrace();
    }
    return null;
  }
  
  
  /**
   * 
   * @return the fatty acid that remain after prediction and cutoff removal
   */
  public Hashtable<String,FattyAcidVO> getAllowedFAs(){
    Hashtable<String,FattyAcidVO> possible = new Hashtable<String,FattyAcidVO>();
    for (FattyAcidVO fa : sortedFAs_.keySet())
      possible.put(fa.getChainId(), fa);
    return possible;
  }
  
  /**
   * 
   * @return the final result; i.e. the relative intensity split according to the prediction - a kind of percentual value
   */
  public Hashtable<String,Double> getRelativeIntensities(){
    return this.relativeIntensities_;
  }
}
