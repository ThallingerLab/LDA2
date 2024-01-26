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

package at.tugraz.genome.lda.msn.vos;

import java.util.Hashtable;

import at.tugraz.genome.lda.utils.LevenbergMarquardtOptimizer;
import at.tugraz.genome.lda.utils.RangeInteger;

/**
 * value object holding information about the current model, the counter model, and the proposed retention time ranges
 * @author Juergen Hartler
 *
 */
public class RtPredictVO
{
  /** the fitted regular LM-model */
  private LevenbergMarquardtOptimizer model_;
  /** the fitted counter LM-model of the hits discarded by MSn*/
  private LevenbergMarquardtOptimizer counterModel_;
  /** the proposed range for the carbon atoms */
  private RangeInteger cAtomsRange_;
  /** the proposed ranges for the double bonds */
  private Hashtable<Integer,RangeInteger> dbsRanges_;
  
  /**
   * constructor requiring all stored objects
   * @param model the fitted regular LM-model
   * @param counterModel fitted counter LM-model of the hits discarded by MSn
   * @param cAtomsRange proposed range for the carbon atoms
   * @param dbsRanges proposed ranges for the double bonds
   */
  public RtPredictVO (LevenbergMarquardtOptimizer model, LevenbergMarquardtOptimizer counterModel, RangeInteger cAtomsRange, Hashtable<Integer,RangeInteger> dbsRanges){
    model_ = model;
    counterModel_ = counterModel;
    cAtomsRange_ = cAtomsRange;
    dbsRanges_ = dbsRanges;
  }
  
  
  /**
   * 
   * @return he fitted regular LM-model
   */
  public LevenbergMarquardtOptimizer getModel()
  {
    return model_;
  }
  /**
   * 
   * @return fitted counter LM-model of the hits discarded by MSn
   */
  public LevenbergMarquardtOptimizer getCounterModel()
  {
    return counterModel_;
  }
  /**
   * 
   * @return proposed range for the carbon atoms
   */
  public RangeInteger getcAtomsRange()
  {
    return cAtomsRange_;
  }
  /**
   * 
   * @return proposed ranges for the double bonds
   */
  public Hashtable<Integer,RangeInteger> getDbsRanges()
  {
    return dbsRanges_;
  }
  
  
}
