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

package at.tugraz.genome.lda.listeners;

import java.io.IOException;

import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.HydroxylationEncodingException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.exception.NoRuleException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.swing.RuleDefinitionInterface;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.quantification.CgException;

/**
 * Listener check the added chain fragment while the user is entering
 * Also paints a new spectra
 * @author Andreas Ziegl
 *
 */
public class AddChainFragmentDocumentListener implements DocumentListener 
{
    private RuleDefinitionInterface rDI_;

    public AddChainFragmentDocumentListener(RuleDefinitionInterface rDI) 
    {
        this.rDI_ = rDI;
    }
    
    public void insertUpdate(DocumentEvent e) 
    {
    	try 
    	{
    		rDI_.checkToAddChainFragment();
    	}
    	catch (RulesException e1) {}
    	catch (IOException e1) {}
    	catch (SpectrummillParserException e1) {}
    	catch (CgException e1) {}
    	catch (NoRuleException e1) {}
    catch (HydroxylationEncodingException | ChemicalFormulaException | LipidCombinameEncodingException e1) {} 
    }
    public void removeUpdate(DocumentEvent e) 
    {
    	try 
    	{
    		rDI_.checkToAddChainFragment();
    	}
    	catch (RulesException e1) {}
    	catch (IOException e1) {}
    	catch (SpectrummillParserException e1) {}
    	catch (CgException e1) {}
    	catch (NoRuleException e1) {} 
    catch (HydroxylationEncodingException | ChemicalFormulaException | LipidCombinameEncodingException e1) {}
    }
    public void changedUpdate(DocumentEvent e) 
    {
    	try 
    	{
    		rDI_.checkToAddChainFragment();
    	}
    	catch (RulesException e1) {}
    	catch (IOException e1) {}
    	catch (SpectrummillParserException e1) {}
    	catch (CgException e1) {}
    	catch (NoRuleException e1) {} 
    catch (HydroxylationEncodingException | ChemicalFormulaException | LipidCombinameEncodingException e1) {}
    }
}
