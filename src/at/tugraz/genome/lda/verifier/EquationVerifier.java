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

package at.tugraz.genome.lda.verifier;

import java.awt.Color;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;

import javax.swing.JOptionPane;
import javax.swing.JTextField;

import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.msn.parser.FragRuleParser;
import at.tugraz.genome.lda.msn.vos.FragmentRuleVO;
import at.tugraz.genome.lda.msn.vos.IntensityRuleVO;
import at.tugraz.genome.lda.swing.RuleDefinitionInterface;

/**
 * Checks if an Equation is correct via the textfield
 * @author Andreas Ziegl
 *
 */
public class EquationVerifier extends JTextField  implements FocusListener
{
  private static final long serialVersionUID = 1L;
  
  /** The textfield to return  */  
  private JTextField myTextField_;
  
  /** The instance of the RDI */
  private RuleDefinitionInterface rdiObject_;
  
  /** The equation string */
  private String equRuleEquationString_;  
  
  /** type 1... head, 2... chain, 3... position */
  int type_;
  
  /** the position in the vector */
  int position_;
  
  /**
   * 
   * @param rdiObject the rule definition interface object
   * @param type 1... head, 2... chain, 3... position
   * @param position the position in the vector
   */  
  public EquationVerifier(RuleDefinitionInterface rdiObject, int type, int position)
  {
    myTextField_ = new JTextField(20);  
    myTextField_.addFocusListener(this); 
    this.rdiObject_ = rdiObject;
    this.type_ = type;
    this.position_ = position;
  }

  /**
   * The aktive textfield to return
   * @return
   */
  public JTextField getTextField()
  {    
    return myTextField_;
  }  
 
  /**
   * Checks if the equation is correct
   */
  public void checkEquation()
  {     
  	IntensityRuleVO ruleVO = null;    
	  try 
	  {
	    if(type_ == 1)
	    {	 
	    	if(rdiObject_.getHeadEquations()[position_].getText() != null && rdiObject_.getHeadEquations()[position_].getText().length()>0)
	     	{ 
	    		this.equRuleEquationString_ = rdiObject_.getHeadEquations()[position_].getText();
	      }
	      if(!(equRuleEquationString_.equals("")))
    	  {
      		ruleVO = FragRuleParser.extractIntensityVOFromEquation(equRuleEquationString_, 0, FragRuleParser.HEAD_SECTION, FragmentRuleVO.getStringKeyHash(rdiObject_.getHeadFragmentRules()), FragmentRuleVO.getStringKeyHash(rdiObject_.getChainFragmentRules()),new Integer(rdiObject_.getAmountOfChains()));
      		ruleVO.setMandatory(rdiObject_.getHeadEquationtMandatories()[position_]);
      		rdiObject_.addEquationToIntensityRules(ruleVO, type_, position_);	
     	  }
	     	else
	     	{	    	    
	     		printStandardEmptyError();     
	     	}
      }
      if(type_ == 2)
      {
        if(rdiObject_.getChainEquations()[position_].getText() != null && rdiObject_.getChainEquations()[position_].getText().length()>0)
        { 
        	this.equRuleEquationString_ = rdiObject_.getChainEquations()[position_].getText();
	      }
	      if(!(equRuleEquationString_.equals("")))
	   	  {
		       ruleVO = FragRuleParser.extractIntensityVOFromEquation(equRuleEquationString_, 0, FragRuleParser.CHAINS_SECTION, FragmentRuleVO.getStringKeyHash(rdiObject_.getHeadFragmentRules()), FragmentRuleVO.getStringKeyHash(rdiObject_.getChainFragmentRules()),new Integer(rdiObject_.getAmountOfChains()));
		       ruleVO.setMandatory(rdiObject_.getChainEquationtMandatories()[position_]);
	         rdiObject_.addEquationToIntensityRules(ruleVO, type_, position_);
	   	  }
	      else
	      {	    	    
	      	printStandardEmptyError();    
	     	}
	    }
	    if(type_ == 3)
	    {
	    	if(rdiObject_.getPositionEquations()[position_].getText() != null && rdiObject_.getPositionEquations()[position_].getText().length()>0)
        { 
        	this.equRuleEquationString_ = rdiObject_.getPositionEquations()[position_].getText();
        }
        if(!(equRuleEquationString_.equals("")))
  	    {
	        ruleVO = FragRuleParser.extractIntensityVOFromEquation(equRuleEquationString_, 0, FragRuleParser.POSITION_SECTION, FragmentRuleVO.getStringKeyHash(rdiObject_.getHeadFragmentRules()), FragmentRuleVO.getStringKeyHash(rdiObject_.getChainFragmentRules()),new Integer(rdiObject_.getAmountOfChains()));
	        ruleVO.setMandatory(rdiObject_.getPositionEquationtMandatories()[position_]);
	        rdiObject_.addEquationToIntensityRules(ruleVO, type_, position_);
	      }  
	    	else
	    	{	    	    
	    		printStandardEmptyError();
	    	}	
	    }
	  }
	  catch (NumberFormatException | RulesException e1) 
    { 	      	
      JOptionPane.showMessageDialog(rdiObject_.topSplitPane_, e1.getMessage());
      myTextField_.requestFocus();
      myTextField_.setBackground(Color.PINK);  
    }	      
   }
  
  /**
   * Prints the error, requests the focus and set the field colored
   */
  public void printStandardEmptyError()
  {
  	JOptionPane.showMessageDialog(rdiObject_.topSplitPane_, "An empty field cannot be saved!"); 
 	  myTextField_.requestFocus();
 	  myTextField_.setBackground(Color.PINK);  
  }
  
  /**
   * The same as check equation but it does not show if there is no success
   */
  public void checkEquationOnlySuccess()
  { 	
    if(!(equRuleEquationString_.equals("")))
    {
      IntensityRuleVO ruleVO = null;
      try 
      {
        if(type_ == 1)
        {        	
          ruleVO = FragRuleParser.extractIntensityVOFromEquation(equRuleEquationString_, 0, FragRuleParser.HEAD_SECTION, FragmentRuleVO.getStringKeyHash(rdiObject_.getHeadFragmentRules()), FragmentRuleVO.getStringKeyHash(rdiObject_.getChainFragmentRules()),new Integer(rdiObject_.getAmountOfChains()));
          ruleVO.setMandatory(rdiObject_.getHeadEquationtMandatories()[position_]);
          rdiObject_.addEquationToIntensityRules(ruleVO, type_, position_);
        }
        if(type_ == 2)
        {        	
          ruleVO = FragRuleParser.extractIntensityVOFromEquation(equRuleEquationString_, 0, FragRuleParser.HEAD_SECTION, FragmentRuleVO.getStringKeyHash(rdiObject_.getHeadFragmentRules()), FragmentRuleVO.getStringKeyHash(rdiObject_.getChainFragmentRules()),new Integer(rdiObject_.getAmountOfChains()));
          ruleVO.setMandatory(rdiObject_.getChainEquationtMandatories()[position_]);
          rdiObject_.addEquationToIntensityRules(ruleVO, type_, position_);
        }
        if(type_ == 3)
        {        	
          ruleVO = FragRuleParser.extractIntensityVOFromEquation(equRuleEquationString_, 0, FragRuleParser.POSITION_SECTION, FragmentRuleVO.getStringKeyHash(rdiObject_.getHeadFragmentRules()), FragmentRuleVO.getStringKeyHash(rdiObject_.getChainFragmentRules()),new Integer(rdiObject_.getAmountOfChains()));
          ruleVO.setMandatory(rdiObject_.getPositionEquationtMandatories()[position_]);
          rdiObject_.addEquationToIntensityRules(ruleVO, type_, position_);
        }
        myTextField_.setBackground(Color.WHITE); 
        rdiObject_.paintNewSpectra(false);
      }
      catch (NumberFormatException | RulesException e1) 
      {
      }      
    }
  }


  /**
   * If the focus is gained on the textfield nothing happens yet
   */
  public void focusGained(FocusEvent e)
  { 
  }

  /**
   * If the focus is lost on the textfield the equation is checked
   */
  public void focusLost(FocusEvent e)
  {
      checkEquation();
  }
  
   /**
    * Initialize the equation string for the first time
    * @param equRuleEquationString
    */
  public void refreshFirstTime()
  {
  	if(type_ == 1)
  	{
  		if(rdiObject_.getHeadEquations()[position_].getText() != null)
  			this.equRuleEquationString_ = rdiObject_.getHeadEquations()[position_].getText();
  	}
  	if(type_ == 2)
  	{
  		if(rdiObject_.getChainEquations()[position_].getText() != null)
  			this.equRuleEquationString_ = rdiObject_.getChainEquations()[position_].getText();
  	}
  	if(type_ == 3)
  	{
  		if(rdiObject_.getPositionEquations()[position_].getText() != null)
  			this.equRuleEquationString_ = rdiObject_.getPositionEquations()[position_].getText();
  	}
  }
  
  /**
   * Refreshes the equation string   
   */
  public void refresh()
  {
  	if(type_ == 1)
  	{
  		if(rdiObject_.getHeadEquations()[position_].getText() != null)
  			this.equRuleEquationString_ = rdiObject_.getHeadEquations()[position_].getText();
  	}
  	if(type_ == 2)
  	{
  		if(rdiObject_.getChainEquations()[position_].getText() != null)
  			this.equRuleEquationString_ = rdiObject_.getChainEquations()[position_].getText();
  	}
  	if(type_ == 3)
  	{
  		if(rdiObject_.getPositionEquations()[position_].getText() != null)
  			this.equRuleEquationString_ = rdiObject_.getPositionEquations()[position_].getText();
  	}
  	
    checkEquationOnlySuccess();   
  }
  
}