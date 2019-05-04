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

package at.tugraz.genome.lda.msn.vos;

import java.util.Hashtable;
import java.util.Vector;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.ShortStringVO;
import at.tugraz.genome.maspectras.parser.spectrummill.ElementConfigParser;
import at.tugraz.genome.maspectras.parser.spectrummill.vos.SmChemicalElementVO;

/**
 * Value object containing all necessary information for a fragment rule
 * @author Juergen Hartler
 *
 */
public class FragmentRuleVO
{
  // identifiers for precursor and chain fragment in the rule file
  public final static String PRECURSOR_NAME = "$PRECURSOR";
  public final static String CHAIN_NAME = "$CHAIN";
  public final static String ALKYL_CHAIN_NAME = "$ALKYLCHAIN";
  public final static String ALKENYL_CHAIN_NAME = "$ALKENYLCHAIN";
  public final static String LCB_NAME = "$LCB";
  
  
  // if a fragment is present, has to be added or subtracted
  private final static int NO_FRAGMENT = 0;
  private final static int ADD_FRAGMENT = 1;
  private final static int MINUS_FRAGMENT = -1;
  

  /** this fragment does not need to be necessarily present*/
  public final static short MANDATORY_UNDEFINED = -1;
  /** this fragment does not need to be necessarily present*/
  public final static short MANDATORY_FALSE = 0;
  /** this fragment must be present, otherwise this hit/structure is discarded*/
  public final static short MANDATORY_TRUE = 1;
  /** this fragment originates from another lipid subclass/adduct*/
  public final static short MANDATORY_OTHER = 2;
  /** this fragment is used for quantitation in PRM data and must be present*/
  public final static short MANDATORY_QUANT = 3;
  /** this fragment is must be present to verify a lipid class - valid for chain fragments only*/
  public final static short MANDATORY_CLASS = 4;
  
  // name of fragment
  private String name_;
  // charge in which the fragment is observed
  private int charge_;
  // MS Level in which the fragment may be found (e.g. MS2, MS3, etc)
  private int msLevel_;
  // must the fragment be present - different options possible according to the MANDATORY_... specifications in this class
  private short mandatory_;
  // is the fragment a neutral loss or adduct to a precursor
  private boolean containsPrecursor_;
  // stati if the chain is present, has to be added or subtracted
  private int chainAction_;
  /** chain type (normal,alkyl,alkenyl)*/  
  private short chainType_;
  // chemical composition of the fragment
  private Hashtable<String,Integer> elementAmounts_;
  //Formula
  private String formula_;
  /** how many hydroxylations must be present for the detection of this fragment; key: number of hydroxylations; value: mandatory - should be null in case of no OH restrictions*/ 
  private RuleHydroxyRequirementSet allowedOHs_;
  /** how many hydroxylations must be present for the detection of this fragment, for a partnering chain; key: number of hydroxylations; value: mandatory - should be null in case of no OH restrictions*/ 
  private RuleHydroxyRequirementSet combiOHs_;

  // the details of the chemical element
  private Hashtable<String,SmChemicalElementVO> elementDetails_;  

  /** the parts of the fragment formula that contains self defined fragments*/
  private Vector<String> selfDefinedParts_;
  
  
  /**
   * @deprecated
   */
  public FragmentRuleVO(String name, String formula, int charge, int msLevel,
      short mandatory, Hashtable<String,FragmentRuleVO> headFragments, 
      Hashtable<String,FragmentRuleVO> chainFragments, ElementConfigParser elementParser) throws RulesException
  {
    this(name,formula,charge,msLevel,mandatory,null,null,headFragments,chainFragments,elementParser);
  }
  
  /**
   * all the information for the VO has to be provided in the constructor
   * @param name name of the fragment
   * @param formula the rule for the fragment
   * @param charge charge state in which the fragment is observed
   * @param msLevel MS Level in which the fragment may be found (e.g. MS2, MS3, etc.)
   * @param mandatory must the fragment be present - different options possible according to the MANDATORY_... specifications in this class
   * @param allowedOHs how many hydroxylations must be present for the detection of this fragment; key: number of hydroxylations; value: mandatory - should be null in case of no OH restrictions
   * @param combiOHs how many hydroxylations must be present for the detection of this fragment; key: number of hydroxylations; value: mandatory - should be null in case of no OH restrictions
   * @param headFragments the previously parsed head fragments (the fragment may be a derivative of a previously parsed head fragment)
   * @param chainFragments the previously parsed chain fragments (the fragment may be a derivative of a previously parsed chain fragment)
   * @param elementParser elementParser for evaluating the chemical formula
   * @throws RulesException specifies in detail which rule has been infringed
   */
  public FragmentRuleVO(String name, String formula, int charge, int msLevel,
      short mandatory, RuleHydroxyRequirementSet allowedOHs, RuleHydroxyRequirementSet combiOHs,
      Hashtable<String,FragmentRuleVO> headFragments, Hashtable<String,FragmentRuleVO> chainFragments,
      ElementConfigParser elementParser) throws RulesException
  {
    super();
    this.name_ = name;
    this.charge_ = charge;
    this.msLevel_ = msLevel;
    this.mandatory_ = mandatory;
    this.formula_ = formula;
    this.chainType_ = LipidomicsConstants.CHAIN_TYPE_NO_CHAIN;
    this.allowedOHs_ = allowedOHs;
    this.combiOHs_ = combiOHs;
    if (charge<1) throw new RulesException("The charge state of an analyte must be greater or equal 1");
    this.categorizeFormula(formula, headFragments, chainFragments, elementParser);
  }

  /**
   * 
   * @return formula
   */
  public String getFormula()
  {
    return formula_;
  }
  
  
  /**
   * 
   * @return name of the fragment
   */
  public String getName()
  {
    return name_;
  }

  /**
   * 
   * @return charge state in which the fragment is observed
   */
  public int getCharge()
  {
    return charge_;
  }

  /**
   * 
   * @return msLevel MS Level in which the fragment may be found (e.g. MS2, MS3, etc.)
   */
  public int getMsLevel()
  {
    return msLevel_;
  }

  /**
   * @param ohNumber the number of hydroxylation sites
   * @return must the fragment be present - different options possible according to the MANDATORY_... specifications in this class
   */
  public short isMandatory(short ohNumber)
  {
    if (ohNumber==LipidomicsConstants.EXCEL_NO_OH_INFO || this.allowedOHs_==null)  
      return mandatory_;
    return (this.allowedOHs_.getEntry(ohNumber).get(0).getMandatory());
  }
  
  /**
   * 
   * @return true when a combiOH is set
   */
  public boolean hasCombiOhRequirements() {
    return this.combiOHs_!=null;
  }
  
  /**
   * @param chainType the type of the partnering chain
   * @param ohNumber the number of hydroxylation sites
   * @return must the fragment be present - different options possible according to the MANDATORY_... specifications in this class
   */
  public short isMandatoryInCombi(short chainType, short ohNumber)
  {
    short mand = FragmentRuleVO.MANDATORY_UNDEFINED;
    if (this.combiOHs_!=null && this.combiOHs_.hasEntry(ohNumber))
      mand = combiOHs_.getMandatory(chainType, ohNumber);
    return mand;
  }
  
  

  public short isMandatory()
  {
    return mandatory_;
  }
  
  
  /**
   * categorizes the composition of the fragment formula (precursor present, chain present or subtracted, other losses or adducts)
   * @param formula the rule for the fragment
   * @param headFragments the previously parsed head fragments (the fragment may be a derivative of a previously parsed head fragment)
   * @param chainFragments the previously parsed chain fragments (the fragment may be a derivative of a previously parsed chain fragment)
   * @param elementParser elementParser for evaluating the chemical formula
   * @throws RulesException specifies in detail which rule has been infringed
   */
  @SuppressWarnings("unchecked")
  private void categorizeFormula(String formula, Hashtable<String,FragmentRuleVO> headFragments, 
      Hashtable<String,FragmentRuleVO> chainFragments, ElementConfigParser elementParser) throws RulesException {
    elementAmounts_ = new Hashtable<String,Integer>();
    elementDetails_ = new Hashtable<String,SmChemicalElementVO>();
    String form = new String(formula);
    Object[] result = containsPrecursor(form);
    form = (String)result[1];
    containsPrecursor_ = (Boolean)result[0];
    result = containsChain(form);
    form = (String)result[1];
    chainAction_ = (Integer)result[0];
    chainType_ = ((Integer)result[2]).shortValue();
    result = containsSelfDefinedFragments(form, getStringKeyHash(headFragments), getStringKeyHash(chainFragments));
    selfDefinedParts_ = (Vector<String>)result[0];
    Vector<Integer> plusOrMinus = (Vector<Integer>)result[1];
    for (int i=0; i!=selfDefinedParts_.size(); i++){
      String selfDefinedName = selfDefinedParts_.get(i); 
      FragmentRuleVO ruleVO = null;
      if (headFragments.containsKey(selfDefinedName)) ruleVO = headFragments.get(selfDefinedName);
      else if (chainFragments.containsKey(selfDefinedName)) ruleVO = chainFragments.get(selfDefinedName);
      int ruleActivity = plusOrMinus.get(i);
      if (ruleVO.containsPrecursor_) containsPrecursor_ = true;
      if (chainFragments.containsKey(selfDefinedName))
        chainType_ = ruleVO.getChainType();
      chainAction_ = chainAction_+ruleActivity*ruleVO.chainAction_;
      for (String element : ruleVO.elementAmounts_.keySet()){
        int elementAmount = 0;
        if (elementAmounts_.containsKey(element)) elementAmount = elementAmounts_.get(element);
        else elementDetails_.put(element, ruleVO.elementDetails_.get(element));
        elementAmount += ruleActivity*ruleVO.elementAmounts_.get(element);
        elementAmounts_.put(element, elementAmount);
      }
    }    
    form = (String)result[2];
    try {
      Hashtable<String,Integer> elements = StaticUtils.categorizeFormula(form);
      for (String element: elements.keySet()){
        if (!elementParser.isElementAvailable(element))
          throw new RulesException("The formula "+formula+" contains the element "+element+" that has not been defined in the "+Settings.getElementConfigPath()+"! Please define the element there, in order to use it!");
        int elementAmount = 0;
        if (elementAmounts_.containsKey(element)) elementAmount = elementAmounts_.get(element);
        else elementDetails_.put(element, elementParser.getElementDetails(element));
        elementAmount += elements.get(element);
        elementAmounts_.put(element, elementAmount);
      }
      if (!elementDetails_.containsKey("H")) elementDetails_.put("H",elementParser.getElementDetails("H"));
      if (!elementDetails_.containsKey("O")) elementDetails_.put("O",elementParser.getElementDetails("O"));
    }
    catch (ChemicalFormulaException e) {
      throw new RulesException("The formula "+formula+" contains fragments that have not been defined before! Fragments have to be defined in previous rows, before they can be used!");
    }
  }
  
  /**
   * checks if the contents of a rule formula is valid
   * @param formula the rule for the fragment
   * @param headFragments the previously parsed head fragments (the fragment may be a derivative of a previously parsed head fragment)
   * @param chainFragments the previously parsed chain fragments (the fragment may be a derivative of a previously parsed chain fragment)
   * @param elementParser elementParser for evaluating the chemical formula
   * @throws RulesException specifies in detail which rule has been infringed
   */
  public static void isFormulaValid (String formula, Hashtable<String,FragmentRuleVO> headFragments, 
      Hashtable<String,FragmentRuleVO> chainFragments, ElementConfigParser elementParser) throws RulesException {
    String form = new String(formula);
    Object[] result = containsPrecursor(form);
    form = (String)result[1];
    result = containsChain(form);
    form = (String)result[1];
    result = containsSelfDefinedFragments(form, getStringKeyHash(headFragments), getStringKeyHash(chainFragments));
    form = (String)result[2];
    try {
      Hashtable<String,Integer> elements = StaticUtils.categorizeFormula(form);
      for (String element: elements.keySet()){
        if (!elementParser.isElementAvailable(element))
          throw new RulesException("The formula "+formula+" contains the element "+element+" that has not been defined in the "+Settings.getElementConfigPath()+"! Please define the element there, in order to use it!");
      }
    }
    catch (ChemicalFormulaException e) {
      throw new RulesException("The formula "+formula+" contains fragments that have not been defined before! Fragments have to be defined in previous columns, before they can be used!");
    }
  }
  
  /**
   * checks if the formula contains $PRECURSOR and removes it from the formula
   * 
   * @param formula
   * @return Object[0] = Boolean (true if $PRECURSOR there); Object[1] = String (original formula reduced by $PRECURSOR)
   * @throws RulesException
   */
  private static Object[] containsPrecursor(String formula) throws RulesException {
    boolean precursorThere = false;
    if (formula.indexOf(PRECURSOR_NAME)!=-1){
      int startIndex = formula.indexOf(PRECURSOR_NAME);
      int stopIndex = startIndex+PRECURSOR_NAME.length();
      if (startIndex!=0){
        startIndex--;
        if (!formula.substring(startIndex,startIndex+1).equalsIgnoreCase("+")) throw new RulesException("Only a \"+\" can be before a "+PRECURSOR_NAME+" in a formula!");
      }
      formula = formula.substring(0,startIndex)+formula.substring(stopIndex,formula.length());
      precursorThere = true;
    }
    Object[] result = new Object[2];
    result[0] = precursorThere;
    result[1] = formula;
    return result;
  }
  
  /**
   * checks if the formula contains $CHAIN and removes it from the formula
   * 
   * @param formula
   * @return Object[0] = Integer (NO_CHAIN/ADD_CHAIN/MINUS_CHAIN; Object[1] = String (original formula reduced by $CHAIN); Object[2] = Integer (ACYL_CHAIN/ALKYL_CHAIN/ALKENYL_CHAIN)
   * @throws RulesException specifies in detail which rule has been infringed
   */
  private static Object[] containsChain(String formula) throws RulesException {
    int chainProcedure = NO_FRAGMENT;
    int chainType = LipidomicsConstants.CHAIN_TYPE_NO_CHAIN;
    if (formula.indexOf(CHAIN_NAME)!=-1 || formula.indexOf(ALKYL_CHAIN_NAME)!=-1 || formula.indexOf(ALKENYL_CHAIN_NAME)!=-1 ||
        formula.indexOf(LCB_NAME)!=-1){
      String chainName = CHAIN_NAME;
      chainType = LipidomicsConstants.CHAIN_TYPE_FA_ACYL;
      if (formula.indexOf(ALKYL_CHAIN_NAME)!=-1){
        chainName = ALKYL_CHAIN_NAME;
        chainType = LipidomicsConstants.CHAIN_TYPE_FA_ALKYL;
      } else if (formula.indexOf(ALKENYL_CHAIN_NAME)!=-1){
        chainName = ALKENYL_CHAIN_NAME;
        chainType = LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL;
      } else if (formula.indexOf(LCB_NAME)!=-1){
        chainName = LCB_NAME;
        chainType = LipidomicsConstants.CHAIN_TYPE_LCB;
      }
      int startIndex = formula.indexOf(chainName);
      int stopIndex = startIndex+chainName.length();
      if (startIndex!=0){
        startIndex--;
        String sign = formula.substring(startIndex,startIndex+1);
        if (sign.equalsIgnoreCase("+")) chainProcedure = ADD_FRAGMENT;
        else if (sign.equalsIgnoreCase("-")) chainProcedure = MINUS_FRAGMENT;
        else throw new RulesException("Only a \"+\" or a \"-\" can be before a "+chainName+" in a formula");
      }else{
        chainProcedure = ADD_FRAGMENT;
      }
      formula = formula.substring(0,startIndex)+formula.substring(stopIndex,formula.length());
    }
    Object[] result = new Object[3];
    result[0] = chainProcedure;
    result[1] = formula;
    result[2] = chainType;
    return result;
  }
  
  /**
   * returns a hash table containing only the keys of a Hashtable<String,FragmentRuleVO> hash
   * @param fragHash the Hashtable<String,FragmentRuleVO>
   * @return hash table containing only the keys of a Hashtable<String,FragmentRuleVO> hash
   */
  public static Hashtable<String,Short> getStringKeyHash(Hashtable<String,FragmentRuleVO> fragHash){
    Hashtable<String,Short> hash = new Hashtable<String,Short>();
    for (String key: fragHash.keySet()) hash.put(key, fragHash.get(key).getChainType());
    return hash;
  }
  
  /**
   * checks if the formula contains any previously defined fragments
   * @param formula chemical formula
   * @param headFragments hash table containing names of already defined head fragments
   * @param chainFragments hash table containing names of already defined chain fragments
   * @return Object[0] = Vector<String> containing used VOs; Object[1] Vector<Integer> telling if chain has to be added or subtracted (ADD_CHAIN/MINUS_CHAIN); Object[1] = String (original formula reduced by the found self defined fragments)
   * @throws RulesException specifies in detail which rule has been infringed
   */
  private static Object[] containsSelfDefinedFragments(String formula, Hashtable<String,Short> headFragments, Hashtable<String,Short> chainFragments) throws RulesException{
    Vector<String> selfDefinedParts = new Vector<String>();
    Vector<Integer> plusOrMinus = new Vector<Integer>();
    Hashtable<String,String> allFragments = new Hashtable<String,String>();
    for (String name : headFragments.keySet()) allFragments.put(name, name);
    for (String name : chainFragments.keySet()) allFragments.put(name, name);
    Vector<ShortStringVO> lengthSortedFragmentNames = getLengthSortedFragmentNames(headFragments,chainFragments);
    
    // if I start with the longest names, I do not have the problem that a smaller name is contained in the name of another fragment
    for (ShortStringVO name: lengthSortedFragmentNames){
      Object[] result = checkSelfDefinedFragment(formula,name.getKey());
      Integer fragmentAction  = (Integer)result[0];
      if (fragmentAction==NO_FRAGMENT) continue;
      selfDefinedParts.add(name.getKey());
      plusOrMinus.add(fragmentAction);
      formula = (String)result[1];
    }
    Object[] result = new Object[3];
    result[0] = selfDefinedParts;
    result[1] = plusOrMinus;
    result[2] = formula;
    return result;
  }
  
  /**
   * sorts fragments names in descending order - this is important for rule interpretation:
   * prevents the detection of smaller fragment names that are contained in longer bigger fragment names
   * instead of the real bigger fragment
   * @param headFragments the names of previously parsed head fragments (the fragment may be a derivative of a previously parsed head fragment)
   * @param chainFragments the names previously parsed chain fragments (the fragment may be a derivative of a previously parsed chain fragment)
   * @return fragment names sorted by length
   */
  public static Vector<ShortStringVO> getLengthSortedFragmentNames(Hashtable<String,Short> headFragments, Hashtable<String,Short> chainFragments){
    return getLengthSortedFragmentNames(headFragments,chainFragments,null);
  }
  
  public static Vector<ShortStringVO> getLengthSortedFragmentNames(Hashtable<String,Short> headFragments, Hashtable<String,Short> chainFragments, Hashtable<String,Short> missed){
    Vector<ShortStringVO> lengthSortedFragmentNames = new Vector<ShortStringVO>();
    for (String name : headFragments.keySet()){
      boolean fragmentAdded = false;
      for (int i=0; i!=lengthSortedFragmentNames.size(); i++){
        if (name.length()>lengthSortedFragmentNames.get(i).getKey().length()){
          lengthSortedFragmentNames.add(i,new ShortStringVO(name,headFragments.get(name)));
          fragmentAdded = true;
          break;
        }
      }
      if (!fragmentAdded) lengthSortedFragmentNames.add(new ShortStringVO(name,headFragments.get(name)));
    }
    for (String name : chainFragments.keySet()){
      boolean fragmentAdded = false;
      for (int i=0; i!=lengthSortedFragmentNames.size(); i++){
        if (name.length()>lengthSortedFragmentNames.get(i).getKey().length()){
          lengthSortedFragmentNames.add(i,new ShortStringVO(name,chainFragments.get(name)));
          fragmentAdded = true;
          break;
        }
      }
      if (!fragmentAdded) lengthSortedFragmentNames.add(new ShortStringVO(name,chainFragments.get(name)));
    }
    if (missed!=null){
      for (String name : missed.keySet()){
        String nm = new String(name);
        if (nm.indexOf("[")!=-1) nm = nm.substring(0,nm.indexOf("["));
        boolean fragmentAdded = false;
        for (int i=0; i!=lengthSortedFragmentNames.size(); i++){
          if (nm.length()>lengthSortedFragmentNames.get(i).getKey().length()){
            lengthSortedFragmentNames.add(i,new ShortStringVO(nm,missed.get(name)));
            fragmentAdded = true;
            break;
          }
        }
        if (!fragmentAdded) lengthSortedFragmentNames.add(new ShortStringVO(nm,missed.get(name)));
      }      
    }
    return lengthSortedFragmentNames;
  }
  
  /**
   * 
   * @param formula formula
   * @param fragment fragment name
   * @return Object[0] = Integer (NO_CHAIN/ADD_CHAIN/MINUS_CHAIN; Object[1] = String (original formula reduced by the found self defined fragments)
   * @throws RulesException specifies in detail which rule has been infringed
   */
  private static Object[] checkSelfDefinedFragment(String formula, String fragment) throws RulesException{
    int fragmentAction = NO_FRAGMENT;
    if (formula.indexOf(fragment)!=-1){
      int startIndex = formula.indexOf(fragment);
      int stopIndex = startIndex+fragment.length();
      if (startIndex!=0){
        startIndex--;
        String sign = formula.substring(startIndex,startIndex+1);
        if (sign.equalsIgnoreCase("+")) fragmentAction = ADD_FRAGMENT;
        else if (sign.equalsIgnoreCase("-")) fragmentAction = MINUS_FRAGMENT;
        else throw new RulesException("Only a \"+\" or a \"-\" can be before the self-defined fragment \""+fragment+"\" in a formula!");
      }else{
        fragmentAction = ADD_FRAGMENT;
      }
      formula = formula.substring(0,startIndex)+formula.substring(stopIndex,formula.length());
    }
    Object[] result = new Object[2];
    result[0] = fragmentAction;
    result[1] = formula;
    return result;
  }
  
  /**
   * @return String representative of the VO - for debugging purposes
   */
  public String toString(){
    String toPrint =   "Name: "+name_;
    toPrint += "; Charge: "+charge_;
    toPrint += "; MS: "+msLevel_;
    toPrint += "; Mand: "+mandatory_;
    toPrint += "; Formula: ";
    if (containsPrecursor_) toPrint += "+"+PRECURSOR_NAME+" ";
    if (chainAction_==ADD_FRAGMENT) toPrint += "+";
    else if (chainAction_==MINUS_FRAGMENT) toPrint += "-";
    if (chainAction_!=NO_FRAGMENT){
      if (chainType_ == LipidomicsConstants.CHAIN_TYPE_FA_ACYL)
        toPrint += CHAIN_NAME+" ";
      else if (chainType_ == LipidomicsConstants.CHAIN_TYPE_FA_ALKYL)
        toPrint += ALKYL_CHAIN_NAME+" ";
      else if (chainType_ == LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL)
        toPrint += ALKENYL_CHAIN_NAME+" ";
      else if (chainType_ == LipidomicsConstants.CHAIN_TYPE_LCB)
        toPrint += LCB_NAME+" ";
    }
    for (String element : elementAmounts_.keySet()){
      if (elementAmounts_.get(element)>0) toPrint += "+";
      else if (elementAmounts_.get(element)<0) toPrint += "-";
      if ((elementAmounts_.get(element)!=0)){
        toPrint += element+Math.abs(elementAmounts_.get(element))+" ";
      }
    }
    if (toPrint.length()>1) toPrint = toPrint.substring(0,toPrint.length()-1);
    return toPrint;
  }

  
  /**
   * calculates the chemical formula of a fragment (by the rules and known precursor and, if necessary, known chain)
   * @param precursorFormula chemical formula of the precursor
   * @param precursorMass m/z of the precursor
   * @param vo
   * @param charge charge state in which the fragment shall be observed
   * @return object vector containing: get(0): chemical formula; get(1) m/z value as double
   * @throws RulesException 
   */
  public Vector<Object> getFormulaAndMass(String precursorFormula, double precursorMass, FattyAcidVO vo, int charge) throws RulesException{
    FattyAcidVO fa = vo;
    return getFormulaAndMass(precursorFormula, precursorMass, fa.getFormula(), fa.getMass(), charge); 
  }

  
  
  /**
   * calculates the chemical formula of a fragment (by the rules and known precursor and, if necessary, known chain)
   * @param precursorFormula chemical formula of the precursor
   * @param precursorMass m/z of the precursor
   * @param chainFormula chemical formula of the analyte chain object
   * @param chainMass mass (not m/z) of the analyte chain object
   * @param charge charge state in which the fragment shall be observed
   * @return object vector containing: get(0): chemical formula; get(1) m/z value as double
   * @throws RulesException 
   */
  public Vector<Object> getFormulaAndMass(String precursorFormula, double precursorMass, String chainFormula, double chainMass, int charge) throws RulesException{   
    Vector<Object> formulaAndMass = new Vector<Object>();
    String formula = "";
    Double mass = 0d;
    Hashtable<String,Integer> formulaAmounts = new Hashtable<String,Integer>();
    if (containsPrecursor_){
      mass = precursorMass;
      try{
        formulaAmounts = StaticUtils.categorizeFormula(precursorFormula);
      } catch (ChemicalFormulaException e) {
        throw new RulesException("The formula "+precursorFormula+" contains fragments that have not been defined before! Fragments have to be defined in previous columns, before they can be used!");
      }
    }
    if (chainAction_!=NO_FRAGMENT){
      Hashtable<String,Integer> chainAmounts = null;
      try{
        chainAmounts = StaticUtils.categorizeFormula(chainFormula);
      } catch (ChemicalFormulaException e) {
        throw new RulesException("The formula "+chainFormula+" contains fragments that have not been defined before! Fragments have to be defined in previous columns, before they can be used!");
      }
      double realChainMass = chainMass;
      if (chainAction_==ADD_FRAGMENT){
        mass += realChainMass/((double)charge);
      }else if (chainAction_==MINUS_FRAGMENT){
        mass -= realChainMass/((double)charge);
      }
      for (String element:chainAmounts.keySet()){
        int amount = 0;
        if (formulaAmounts.containsKey(element)) amount = formulaAmounts.get(element);
        if (chainAction_==ADD_FRAGMENT) amount += chainAmounts.get(element);
        else if (chainAction_==MINUS_FRAGMENT) amount -= chainAmounts.get(element);
        formulaAmounts.put(element, amount);
      }
    }
    for (String element : elementAmounts_.keySet()){
      int amount = 0;
      if (formulaAmounts.containsKey(element)) amount = formulaAmounts.get(element);
      int elementAmount = elementAmounts_.get(element);
      amount += elementAmount;
      formulaAmounts.put(element, amount);
      mass += (elementAmount*elementDetails_.get(element).getMonoMass())/((double)charge);
    }
    for (String element : formulaAmounts.keySet()){
      if (formulaAmounts.get(element)>0) formula += "+";
      else if (formulaAmounts.get(element)<0) formula += "-";
      if ((formulaAmounts.get(element)!=0)){
        formula += element+Math.abs(formulaAmounts.get(element))+" ";
      }
    }
    formulaAndMass.add(formula);
    formulaAndMass.add(mass);
    return formulaAndMass;
  }

  /**
   * 
   * @return the type of the chain (ACYL_CHAIN/ALKYL_CHAIN/ALKENYL_CHAIN/LCB)
   */
  public short getChainType()
  {
    return chainType_;
  }
  
  /**
   * checks if this rule VO contains valid fragments only (if one was deleted, the depending rule)
   * has to be deleted as well
   * @param allowedRules the fragment rule names that are valid
   * @return true if it contains only fragments that are valid
   */
  public boolean containsOnlyAllowedFragments(Hashtable<String,FragmentRuleVO> allowedRules){
    boolean containsOnlyAllowed = true;
    for (String fragment : selfDefinedParts_){
      if (!allowedRules.containsKey(fragment)){
        containsOnlyAllowed = false;
        break;
      }
    }
    return containsOnlyAllowed;
  }
  
  
  /**
   * checks whether this fragment is applicable for this OH configuration
   * @param ohNumber the number of hydroxylation sites
   * @return true when the fragment is applicable
   */
  public boolean hydroxylationValid(short ohNumber) {
    if (ohNumber==LipidomicsConstants.EXCEL_NO_OH_INFO || this.allowedOHs_==null)
      return true;
    return (this.allowedOHs_.hasEntry(ohNumber));
  }
  
}
