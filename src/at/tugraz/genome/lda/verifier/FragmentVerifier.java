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



import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JTextField;

import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.msn.vos.FragmentRuleVO;
import at.tugraz.genome.lda.swing.RuleDefinitionInterface;
import at.tugraz.genome.lda.swing.rdi.FragmentFormulaTextField;
import at.tugraz.genome.lda.swing.rdi.FragmentNameTextField;


/**
 * Checks if the input for fragments is correct via the textfields
 * @author Andreas Ziegl
 *
 */
public class FragmentVerifier
{
	
/** The text field to return */
  private JTextField myTextField_;
  
  /** The instance of the RDI */
  private RuleDefinitionInterface rdiObject_;
  
  /** Name of the fragment */
  private String fragRuleNameString_;
  
  /** Formula of the fragment */
  private String fragRuleFormulaString_;
  
  /** Charge of the fragment */
  private String fragRuleChargeString_;
  
  /** MS Level of the fragment */
  private String fragRuleMsLevelString_;
  
  /** The mandatory of the fragment */
  private boolean mandatory_;
  
  /** 1... formula, 2... charge, 3... ms level */
  int type_;
  
  /** Shows if it is a head or chain fragment */
  int headOrChain_;
  
  /** The position of the element in the table */
  int position_;
   
  public static int TYPE_NAME = 0;
  public static int TYPE_FORMULA = 1;
  public static int TYPE_CHARGE = 2;
  public static int TYPE_MSLEVEL = 3;
  
  
  /**
   * 
   * @param rdiObject the rule definition interface object
   * @param type 1... formula, 2... charge, 3... ms level
   * @param headOrChain 1... head 2... chain
   */  
  public FragmentVerifier(RuleDefinitionInterface rdiObject, int type, int headOrChain)
  { 
    this.type_ = type;
    this.createTextFieldVariable();
    this.rdiObject_ = rdiObject;
    this.headOrChain_ = headOrChain;
  }
  
  /**
   * Creates the textfile for the verification
   */
  public void createTextFieldVariable()
  {
    if (type_ ==TYPE_NAME)
      myTextField_ = new FragmentNameTextField();
    else if (type_==TYPE_FORMULA)
      myTextField_ = new FragmentFormulaTextField();
    else
  	    myTextField_ = new JTextField(20);
  }

  /**
   * Returns the active textfield
   * @return the JTextField
   */
  public JTextField getTextField()
  {    
    return this.myTextField_;
  }  
  
  /**
   * Checks the fragment and saves the txt in the cache folder to refresh the container for the spectra
   */
  public void checkFragment()
  {
  	if(rdiObject_.errorMessageNeeded_)
  	{
			rdiObject_.errorMessageNeeded_ = false;			
			if(headOrChain_ == 1)
			{
				this.fragRuleNameString_ = rdiObject_.getHeadFragmentVector().elementAt(position_).getName(); 		
				this.fragRuleFormulaString_ = rdiObject_.getHeadRuleFormulas()[position_].getText();		
			  this.fragRuleChargeString_ = rdiObject_.getHeadRuleCharges()[position_].getText();		
			  this.fragRuleMsLevelString_ = rdiObject_.getHeadRuleMSLevels()[position_].getText();			  
			  this.mandatory_ = rdiObject_.getHeadFragmentMandatories()[position_];			  
			}
			
			if(headOrChain_ == 2)
			{
				this.fragRuleNameString_ = rdiObject_.getChainFragmentVector().elementAt(position_).getName(); 		
				this.fragRuleFormulaString_ = rdiObject_.getChainRuleFormulas()[position_].getText();		
			  this.fragRuleChargeString_ = rdiObject_.getChainRuleCharges()[position_].getText();		
			  this.fragRuleMsLevelString_ = rdiObject_.getChainRuleMSLevels()[position_].getText();	
			  this.mandatory_ = rdiObject_.getChainFragmentMandatories()[position_];	
			}
			checkMsLevel(true);
				if(!(fragRuleFormulaString_.equals("")))
				{
				  try 
				  {
				    FragmentRuleVO ruleVO = null;
					  if(type_ == 1)
					  {
					    ruleVO = rdiObject_.checkFragment(fragRuleNameString_, fragRuleFormulaString_, 0, 0, mandatory_, 1, headOrChain_);  
					  }
			      
					  
					  if(!(fragRuleChargeString_.equals("")))
					  {
						  if(type_ == 2)
						  {
						    ruleVO = rdiObject_.checkFragment(fragRuleNameString_, fragRuleFormulaString_, Integer.parseInt(fragRuleChargeString_), 0, mandatory_, 2, headOrChain_); 	  
						  }
					  
			      
						  if(!(fragRuleMsLevelString_.equals("")))
						  {
							  if(type_ == 3)
							  {
							  	checkMsLevel(true);
							  ruleVO = rdiObject_.checkFragment(fragRuleNameString_, fragRuleFormulaString_, Integer.parseInt(fragRuleChargeString_), Integer.parseInt(fragRuleMsLevelString_), mandatory_, 3, headOrChain_);    	  
							  }   
						  }
					  }
					    if (ruleVO!=null) rdiObject_.updateRuleVO(ruleVO,headOrChain_);
				  }
			    catch (RulesException e1) 
			    {
	      myTextField_.requestFocus();
			      JOptionPane.showMessageDialog(rdiObject_.topSplitPane_, e1.getMessage()+ " The last accepted value remains in the textfile.");			      
			      myTextField_.setBackground(Color.PINK);
			    }               
			    catch (java.lang.NumberFormatException e2)
			    {
			    	myTextField_.setBackground(Color.PINK);
			      JOptionPane.showMessageDialog(rdiObject_.topSplitPane_, "Not a number! The last accepted value remains in the textfile."); 
			      myTextField_.requestFocus();			      
			    }
				
			}
			else
			{
				myTextField_.requestFocus();
			  JOptionPane.showMessageDialog(rdiObject_.topSplitPane_, "An empty field cannot be saved! The last accepted value remains in the textfile.");
			  
			  myTextField_.setBackground(Color.PINK);
			}
  	}   
  }
  
  /**
   * Checks without exception to not interrupt the user during the writing
   */
  public void checkFragmentOnlySuccess()
  {
  	if(!checkMsLevel(false))
  	{
  		rdiObject_.errorMessageNeeded_ = true;
  	}
  	else
  	{
  		myTextField_.setBackground(Color.WHITE);
  	}
  	if(!(fragRuleFormulaString_.equals("")))
		{
		  try 
		  {
		    FragmentRuleVO ruleVO = null;
			  if(type_ == 1)
			  {
			    ruleVO = rdiObject_.checkFragment(fragRuleNameString_, fragRuleFormulaString_, 0, 0, mandatory_, 1, headOrChain_);  
			  }
	      
			  
			  if(!(fragRuleChargeString_.equals("")))
			  {
				  if(type_ == 2)
				  {
				    ruleVO = rdiObject_.checkFragment(fragRuleNameString_, fragRuleFormulaString_, Integer.parseInt(fragRuleChargeString_), 0, mandatory_, 2, headOrChain_); 	  
				  }
			  
	      
				  if(!(fragRuleMsLevelString_.equals("")))
				  {
					  if(type_ == 3)
					  {
					    ruleVO = rdiObject_.checkFragment(fragRuleNameString_, fragRuleFormulaString_, Integer.parseInt(fragRuleChargeString_), Integer.parseInt(fragRuleMsLevelString_), mandatory_, 3, headOrChain_);    	  
					  }   
				  }
				    if (ruleVO!=null) rdiObject_.updateRuleVO(ruleVO,headOrChain_);
			  }
		  }
    catch (RulesException e1) 
    { 
    	rdiObject_.errorMessageNeeded_ = true;
    }               
    catch (java.lang.NumberFormatException e2)
    {   
    	rdiObject_.errorMessageNeeded_ = true;
    }
		}
	  else
	  {
		  rdiObject_.errorMessageNeeded_ = true;
	  }
  }

  /**
   * Check Ms Level
   * 
   */
  public boolean checkMsLevel(boolean errormessage)
  {
  	if(fragRuleMsLevelString_.equals(""))
  	{
  		if(errormessage == true)
  		{
  			myTextField_.requestFocus();
  		  JOptionPane.showMessageDialog(rdiObject_.topSplitPane_, "This ms level is not possible.");
  		  this.fragRuleMsLevelString_ = "2";
  		  myTextField_.setBackground(Color.PINK);
  		}  		
		  return false;
  	}
  	else
  	{
    	if(!(fragRuleMsLevelString_.equals("2")))
    	{
    		if(errormessage == true)
    		{
	    		myTextField_.requestFocus();
	  		  JOptionPane.showMessageDialog(rdiObject_.topSplitPane_, "This ms level is not possible.");
	  		  this.fragRuleMsLevelString_ = "2";
	  		  myTextField_.setBackground(Color.PINK);
    		}
  		  return false;
    	}
    	else
    	{     		
    		return true;
    	}
  	}
  
  }
  
  /**
   * Initializes the variables without refreshing the spectra
   * @param position
   */
  public void refreshFirstTime(int position)
  {		
  	this.position_ = position;
  	if(headOrChain_ == 1)
  	{
	    this.fragRuleNameString_ = rdiObject_.getHeadFragmentVector().elementAt(position).getName();  
	    this.fragRuleFormulaString_ = rdiObject_.getHeadRuleFormulas()[position].getText();   
	    this.fragRuleChargeString_ = rdiObject_.getHeadRuleCharges()[position].getText();   
	    this.fragRuleMsLevelString_ = rdiObject_.getHeadRuleMSLevels()[position].getText();	   
	    this.mandatory_ = rdiObject_.getHeadFragmentMandatories()[position];
  	}
  	if(headOrChain_ == 2)
  	{
	    this.fragRuleNameString_ = rdiObject_.getChainFragmentVector().elementAt(position).getName();  
	    this.fragRuleFormulaString_ = rdiObject_.getChainRuleFormulas()[position].getText();   
	    this.fragRuleChargeString_ = rdiObject_.getChainRuleCharges()[position].getText();   
	    this.fragRuleMsLevelString_ = rdiObject_.getChainRuleMSLevels()[position].getText();
	    this.mandatory_ = rdiObject_.getChainFragmentMandatories()[position];	
  	}
  }  

  
/**
 * If something in the textfield changes this function starts to work
 * @param fragRuleNameString
 * @param fragRuleFormulaString
 * @param fragRuleChargeString
 * @param fragRuleMsLevelString
 */
  public void refresh(int position)
  {		
  	this.position_ = position;
  	if(headOrChain_ == 1)
  	{
			if(rdiObject_.getHeadFragmentVector().elementAt(position).getName() != null && rdiObject_.getHeadFragmentVector().elementAt(position).getName().length()>0)
				this.fragRuleNameString_ = rdiObject_.getHeadFragmentVector().elementAt(position).getName(); 
			
			if(checkFragmentFormulaInput(rdiObject_.getHeadRuleFormulas()[position]))
				this.fragRuleFormulaString_ = rdiObject_.getHeadRuleFormulas()[position].getText(); 
			
			if(rdiObject_.getHeadRuleCharges()[position].getText() != null && rdiObject_.getHeadRuleCharges()[position].getText().length()>0)
		    	this.fragRuleChargeString_ = rdiObject_.getHeadRuleCharges()[position].getText();  
			
			if(rdiObject_.getHeadRuleMSLevels()[position].getText() != null && rdiObject_.getHeadRuleMSLevels()[position].getText().length()>0)
			    	this.fragRuleMsLevelString_ = rdiObject_.getHeadRuleMSLevels()[position].getText();

	    this.mandatory_ = rdiObject_.getHeadFragmentMandatories()[position];		    
  	}
  	if(headOrChain_ == 2)
   	{
			if(rdiObject_.getChainFragmentVector().elementAt(position).getName() != null && rdiObject_.getChainFragmentVector().elementAt(position).getName().length()>0)
				this.fragRuleNameString_ = rdiObject_.getChainFragmentVector().elementAt(position).getName(); 
			
			if(checkFragmentFormulaInput(rdiObject_.getChainRuleFormulas()[position]))
				this.fragRuleFormulaString_ = rdiObject_.getChainRuleFormulas()[position].getText(); 
			
			if(rdiObject_.getChainRuleCharges()[position].getText() != null && rdiObject_.getChainRuleCharges()[position].getText().length()>0)
		    	this.fragRuleChargeString_ = rdiObject_.getChainRuleCharges()[position].getText();  
			
			if(rdiObject_.getChainRuleMSLevels()[position].getText() != null && rdiObject_.getChainRuleMSLevels()[position].getText().length()>0)
		    	this.fragRuleMsLevelString_ = rdiObject_.getChainRuleMSLevels()[position].getText();
			
			this.mandatory_ = rdiObject_.getChainFragmentMandatories()[position];			
  	}
    this.checkFragmentOnlySuccess();  
    if(checkMsLevel(false))
    {
    	rdiObject_.paintNewSpectra(false); 
    }       
  }
  
  /**
   * for checking a JTextField of the type MSnAnalyzer fragment annotation name if it contains any not allowed characters
   * @param field the JTextField
   * @return true if the field is OK
   */
  public static boolean checkFragmentNameInput(JTextField field){
    boolean ok = true;
    if (field==null || field.getText()==null || field.getText().length()==0) return false;
    String text = field.getText().trim();
    if (text.indexOf(" ")!=-1){
      new WarningMessage(new JFrame(), "Error", "A rule name must not contain any empty spaces");
      text.replaceAll(" ", "");
      field.setText(text);
    }
    if (text.length()==0) return false;
    if (text.indexOf("=")!=-1){
      new WarningMessage(new JFrame(), "Error", "A rule name must not contain any \"=\" signs");
      text.replaceAll("=", "");
      field.setText(text);
    }
    if (text.length()==0) return false;
    return ok;
  }
  
  /**
   * for checking a JTextField of the type MSnAnalyzer fragment formula if it contains any not allowed characters
   * @param field the JTextField
   * @return true if the field is OK
   */
  public static boolean checkFragmentFormulaInput(JTextField field){
    boolean ok = true;
    if (field==null || field.getText()==null || field.getText().length()==0) return false;
    String text = field.getText().trim();
    if (text.indexOf(" ")!=-1){
      new WarningMessage(new JFrame(), "Error", "A rule formula must not contain any empty spaces");
      return false;
    }
    if (text.length()==0) return false;
    if (text.indexOf("=")!=-1){
      new WarningMessage(new JFrame(), "Error", "A rule formula must not contain any \"=\" signs");
      return false;
    }
    if (text.length()==0) return false;
    return ok;
  }
  

}