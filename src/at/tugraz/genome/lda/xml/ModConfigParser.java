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

package at.tugraz.genome.lda.xml;

import javax.xml.parsers.DocumentBuilderFactory;  
import javax.xml.parsers.DocumentBuilder;  
import org.w3c.dom.Document;  
import org.w3c.dom.NodeList;


import org.w3c.dom.Node;  
import java.io.File; 
import java.util.Hashtable;

/**
 * 
 * @author Christoph Krettler
 *
 */
public class ModConfigParser {
	
	private Hashtable<String,Hashtable<String,Integer>> modifications_;
	private File modConfigFile_;
	
	
	/**
	  * constructor 
	  * @param modConfigPath path of the modconfig.xml file
	*/
	public ModConfigParser(String modConfigPath_){
		this.modConfigFile_ = new File(modConfigPath_);
	}
	
	/**
	  * parses a given  modconfig.xml file
	  * @param modConfigPath path of the modconfig.xml file
	*/
	public void parse()
	{
		modifications_ = new Hashtable<String,Hashtable<String,Integer>>(); 
		try   
		{  
			DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();  
			DocumentBuilder db = dbf.newDocumentBuilder();  
			Document doc = db.parse(modConfigFile_);  
			doc.getDocumentElement().normalize();  
			
			NodeList nodeList = doc.getElementsByTagName("modification");  
			for (int i = 0; i < nodeList.getLength(); i++)   
			{  
				Node node = nodeList.item(i);
				String symbol = node.getAttributes().getNamedItem("symbol").getNodeValue();
				Hashtable<String,Integer> elements = new Hashtable<String,Integer>();
				
				NodeList subNodeList = node.getChildNodes();  
				for (int i2 = 0; i2 < subNodeList.getLength(); i2++)
				{
					Node subNode = subNodeList.item(i2);
					if(subNode.hasAttributes())
					{
						String element = subNode.getAttributes().getNamedItem("name").getNodeValue();
						int amount = Integer.parseInt(subNode.getTextContent());
						elements.put(element,amount);
					}
					
				}
				modifications_.put(symbol, elements);
			}  
		}   
		catch (Exception e)   
		{  
			e.printStackTrace();  
		}  
	}
	
	/**
	  * returns the change in elements for a given modification
	  * @param modification the modification for which the composition wants to be known
	  * @return elements and their amount (can be negative)  
	*/
	public Hashtable<String,Integer> getElements(String modification) throws Exception
	{
		if(modifications_.containsKey(modification))
		{
			return modifications_.get(modification);
		}
		else {
			throw new Exception("The modification " + modification + " is listed in the fattyAcidChains file but not in the modconfig.xml file! Will continue withouth it...");
		}
		
	}
}
