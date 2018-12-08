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

package at.tugraz.genome.lda.vos.rdi;

import at.tugraz.genome.lda.msn.RulesContainer;

/**
 * VO for storing values of the general settings input tab of the rule definition interface (RDI)
 * @author Juergen Hartler
 *
 */
public class GeneralSettingsVO
{
  /** how many chains does this analyte class have */
  private Integer amountOfChains_;
  /** how many alkyl chains does this analyte class have */
  private Integer amountOfAlkylChains_;
  /** how many alkenyl chains does this analyte class have */
  private Integer amountOfAlkenylChains_;
  /** how many LCBs does this analyte class have */
  private Short amountOfLCBs_;
  /** the name of the fatty acid chain library */
  private String chainLibrary_;
  /** the name of the long chain base library */
  private String lcbLibrary_;
  /**  Java regular expression to extract the number carbon atoms from the analyte name */
  private String carbonAtomsRule_;
  /** Java regular expression to extract the number of double bonds from the analyte name */
  private String doubleBondsRule_;
  /** is an identification which is based on a single chain already valid? */
  private boolean allowSingleChain_;
  /** cutoff value relative to the highest chain combination found */
  private String chainCutoff_;
  /** cutoff value relative to the base peak (in its original format)*/
  private String basePeakCutoff_;
  /** spectrum coverage that has to be fulfilled by the found fragments (originally value including % or permille sign)*/
  private String spectrumCoverage_;
  /** shall for this class a post processing by retention time be executed */
  private boolean rtPostProcessing_;
  /** shall the retention time processing take a potential parallel model into account */
  private boolean rtParallelSeries_;
  /** the maximum retention time which a value may deviate from the hypothetical RT model */
  private Double rtMaxDeviation_;
  /** in which order the identification by MS should be made (ORDER_MS1_FIRST/ORDER_MSN_FIRST/ORDER_MSN_ONLY from RulesContainer)*/
  private int msIdentificationOrder_;
  /** are there allowed positions*/
  private Integer addChainPositions_;
  
  
  /**
   * Constructor setting the VO to the empty default values (should be called only if a new rule set is entered)
   */
  public GeneralSettingsVO(){
    amountOfChains_ = null;
    amountOfAlkylChains_ = null;
    amountOfAlkenylChains_ = null;
    amountOfLCBs_ = null;
    chainLibrary_ = null;
    lcbLibrary_ = null;
    carbonAtomsRule_ = null;
    doubleBondsRule_ = null;
    allowSingleChain_ = false;
    chainCutoff_ = "default";
    basePeakCutoff_ = "0%";
    spectrumCoverage_ = "0%";
    rtPostProcessing_ = true;
    rtParallelSeries_ = false;
    rtMaxDeviation_ = null;
    msIdentificationOrder_ = RulesContainer.ORDER_MS1_FIRST;
  }
  
  /**
   * constructor setting the parameters of the VO correspondingly (for existing fragmentation rule sets)
   * @param amountOfChains how many chains does this analyte class have
   * @param amountOfAlkylChains how many alkyl chains does this analyte class have
   * @param amountOfAlkenylChains how many alkenyl chains does this analyte class have
   * @param amountOfLCBs how many LCBs does this analyte class have
   * @param chainLibrary the name of the fatty acid chain library
   * @param chainLibrary the name of the long chain base library
   * @param carbonAtomsRule Java regular expression to extract the number carbon atoms from the analyte name
   * @param doubleBondsRule Java regular expression to extract the number of double bonds from the analyte name
   * @param allowSingleChain is an identification which is based on a single chain already valid?
   * @param chainCutoff cutoff value relative to the highest chain combination found
   * @param basePeakCutoff cutoff value relative to the base peak (in its original format)
   * @param spectrumCoverage spectrum coverage that has to be fulfilled by the found fragments (originally value including % or permille sign)
   * @param rtPostProcessing  shall for this class a post processing by retention time be executed
   * @param rtParallelSeries shall the retention time processing take a potential parallel model into account
   * @param rtMaxDeviation the maximum retention time which a value may deviate from the hypothetical RT model
   * @param msIdentificationOrder in which order the identification by MS should be made (ORDER_MS1_FIRST/ORDER_MSN_FIRST/ORDER_MSN_ONLY from RulesContainer)
   */
  public GeneralSettingsVO(Integer amountOfChains,
      Integer amountOfAlkylChains, Integer amountOfAlkenylChains, Short amountOfLCBs,
      Integer addChainPositions, String chainLibrary, String lcbLibrary, String carbonAtomsRule,
      String doubleBondsRule, boolean allowSingleChain, String chainCutoff, String basePeakCutoff,
      String spectrumCoverage, boolean rtPostProcessing, boolean rtParallelSeries,
      Double rtMaxDeviation, int msIdentificationOrder)
  {
    this();
    this.amountOfChains_ = amountOfChains;
    this.amountOfAlkylChains_ = amountOfAlkylChains;
    this.amountOfAlkenylChains_ = amountOfAlkenylChains;
    this.amountOfLCBs_ = amountOfLCBs;
    this.addChainPositions_ = addChainPositions;
    this.chainLibrary_ = chainLibrary;
    this.lcbLibrary_ = lcbLibrary;
    this.carbonAtomsRule_ = carbonAtomsRule;
    this.doubleBondsRule_ = doubleBondsRule;
    this.allowSingleChain_ = allowSingleChain;
    this.chainCutoff_ = chainCutoff;
    this.basePeakCutoff_ = basePeakCutoff;
    this.spectrumCoverage_ = spectrumCoverage;
    this.rtPostProcessing_ = rtPostProcessing;
    this.rtParallelSeries_ = rtParallelSeries;
    this.rtMaxDeviation_ = rtMaxDeviation;
    this.msIdentificationOrder_ = msIdentificationOrder;
  }

  /**
   * 
   * @return how many chains does this analyte class have
   */
  public Integer getAmountOfChains()
  {
    return amountOfChains_;
  }

  /**
   * 
   * @return how many alkyl chains does this analyte class have
   */
  public Integer getAmountOfAlkylChains()
  {
    return amountOfAlkylChains_;
  }

  /**
   * 
   * @return how many alkenyl chains does this analyte class have
   */
  public Integer getAmountOfAlkenylChains()
  {
    return amountOfAlkenylChains_;
  }
  
  /**
   * 
   * @return how many LCBs does this analyte class have
   */
  public Short getAmountOfLCBs()
  {
    return amountOfLCBs_;
  }

  /**
   * 
   * @return the name of the fatty acid chain library
   */
  public String getChainLibrary()
  {
    return chainLibrary_;
  }

  /**
   * 
   * @return the name of the long chain base library
   */
  public String getLcbLibrary()
  {
    return lcbLibrary_;
  }

  /**
   * 
   * @return Java regular expression to extract the number carbon atoms from the analyte name
   */
  public String getCarbonAtomsRule()
  {
    return carbonAtomsRule_;
  }

  /**
   * 
   * @return Java regular expression to extract the number of double bonds from the analyte name
   */
  public String getDoubleBondsRule()
  {
    return doubleBondsRule_;
  }

  /**
   * 
   * @return is an identification which is based on a single chain already valid?
   */
  public boolean isAllowSingleChain()
  {
    return allowSingleChain_;
  }

  /**
   * 
   * @return cutoff value relative to the highest chain combination found
   */
  public String getChainCutoff()
  {
    return chainCutoff_;
  }

  /**
   * 
   * @return cutoff value relative to the base peak (in its original format)
   */
  public String getBasePeakCutoff()
  {
    return basePeakCutoff_;
  }

  /**
   * 
   * @return spectrum coverage that has to be fulfilled by the found fragments (originally value including % or permille sign)
   */
  public String getSpectrumCoverage()
  {
    return spectrumCoverage_;
  }

  /**
   * 
   * @return shall for this class a post processing by retention time be executed
   */
  public boolean isRtPostProcessing()
  {
    return rtPostProcessing_;
  }

  /**
   * 
   * @return shall the retention time processing take a potential parallel model into account
   */
  public boolean isRtParallelSeries()
  {
    return rtParallelSeries_;
  }

  /**
   * 
   * @return the maximum retention time which a value may deviate from the hypothetical RT model
   */
  public Double getRtMaxDeviation()
  {
    return rtMaxDeviation_;
  }

  /**
   * 
   * @return in which order the identification by MS should be made (ORDER_MS1_FIRST/ORDER_MSN_FIRST/ORDER_MSN_ONLY from RulesContainer)
   */
  public int getMsIdentificationOrder()
  {
    return msIdentificationOrder_;
  }

  /**
   * 
   * @return are there allowed positions
   */
  public Integer getAddChainPositions()
  {
    return addChainPositions_;
  }

  
}
