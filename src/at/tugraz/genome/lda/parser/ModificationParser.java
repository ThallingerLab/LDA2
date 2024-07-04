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

package at.tugraz.genome.lda.parser;

import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import at.tugraz.genome.lda.Settings;

/**
 * 
 * @author Christoph Krettler
 *
 */
public class ModificationParser {

	Hashtable<String, Integer> modFormula_;
	String modification_;
	
	/**
	  * constructor 
	  * @param modification the modification to parse
	*/
	public ModificationParser(String modification) {
		this.modification_ = modification;
	}
	
	/**
	  * parses the modification and returns elemental change (O, 4OH) --> -2H + 4O 
	  * @param modification the modification to parse
	*/
	public void parse() {
		Hashtable<String, Integer> modFormula = new Hashtable<String, Integer>();
		
		if(modification_.equals(""))
		{
			modFormula_ = modFormula;
			return;
		}
		
    	for(String mod : modification_.split(","))
    	{
    		int multiplier = 1;
    		
    //		if(Character.isDigit(mod.charAt(0)))
    //		{
    			Pattern p = Pattern.compile("(\\d+)");
        	    Matcher m = p.matcher(mod);
        	    if(m.find())
        	    {
        	    	multiplier = Integer.parseInt(m.group(0));
        	    	mod = mod.replace(m.group(0), "");
        	    }         	   
    	//	}
    		try {
				for (String element : Settings.getModParser().getElements(mod).keySet())
				{
					if(!modFormula.containsKey(element))
					{
						modFormula.put(element, multiplier*Settings.getModParser().getElements(mod).get(element));
					}
					else
					{
						modFormula.put(element,modFormula.get(element)+multiplier*Settings.getModParser().getElements(mod).get(element));
					}
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
    	}
    	modFormula_ = modFormula;
	}
	
	/**
	  * Adds modifications to a chemical formula, returns new formula as string 
	  * @param elementalComposition the old elemental composition of a formula 
	  * @return the new elemental composition of a formula (+added modification; as string)
	*/
	public String getNewChemicalFormula(Hashtable<String, Integer> elementalComposition){
		String chemicalFormula = "";
		for (String element : elementalComposition.keySet()){
	          if (chemicalFormula.length()>0) chemicalFormula+=" ";
	          int amount = elementalComposition.get(element);
	          
	          if(modFormula_.containsKey(element))
	          {
	        	  chemicalFormula+=element+String.valueOf(amount+modFormula_.get(element));
	          }
	          else
	          {
	        	  chemicalFormula+=element+String.valueOf(amount);
	          }
	     }
	     for (String element : modFormula_.keySet()){
		     if(!elementalComposition.containsKey(element) ) {
		    	 chemicalFormula+=element+String.valueOf( modFormula_.get(element));
		     }
	     }
	     return chemicalFormula;
	}
	
	/**
	  * Adds modifications to a chemical composition, returns new composition as Hashtable
	  * @param elementalComposition the old elemental composition of a formula 
	  * @return the new elemental composition of a formula (+added modification; as a Hashtable)
	*/
	public Hashtable<String, Integer> getNewChemicalComposition(Hashtable<String, Integer> elementalComposition){
		Hashtable<String, Integer> chemicalComposition = new Hashtable<String, Integer>();
		for (String element : elementalComposition.keySet()){
	          if(modFormula_.containsKey(element))
	          {
	        	  chemicalComposition.put(element, elementalComposition.get(element) + modFormula_.get(element));
	          }
	          else
	          {
	        	  chemicalComposition.put(element, elementalComposition.get(element));
	          }
	     }
	     for (String element : modFormula_.keySet()){
		     if(!elementalComposition.containsKey(element) ) {
		    	 chemicalComposition.put(element, modFormula_.get(element));
		     }
	     }
		return chemicalComposition;
	     
	}
	
	/**
	  * Adds modifications to a mass
	  * @param mass the old mass 
	  * @return the new mass (+added modification)
	*/
	public double getNewMass(double mass) {
		for (String element : modFormula_.keySet())
		{
			mass += modFormula_.get(element)*Settings.getElementParser().getElementDetails(element).getMonoMass();
		}
		return mass;
	}
	
	/**
	  * Returns chemical composition of a given modification (as a string)
	  * @return the chemical composition of a given modification (as a string)
	*/
	public String getModificationComposition() {
		return this.modFormula_.toString();
	}

}
