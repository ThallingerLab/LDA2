package at.tugraz.genome.lda.masslist;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.util.ArrayList;

import javax.swing.JButton;
import javax.swing.JPanel;

import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.vos.AdductVO;

/**
 * 
 * @author Leonida M. Lamp
 * 
 * 
 * We need... a folder with basic values for different classes
 * This would be... oh number, oh range, rt start, rt stop, adductinsensitivertfilter, bestMatchBySpectrumCoverage
 * (Also ox stuff)
 * Adducts (save in separate folder ..name,formula,charge)
 * Formula of headgroup (without chains)
 * 
 * Minimum number of chain carbon atoms
 * Maximum number of chain carbon atoms
 * Minimum number of chain double bonds
 * Maximum number of chain double bonds
 * 
 * General setting: SIL? path to fattyacidchainlist.xlsx, only if SIL desired (to compute possible combinations based on given params).
 * 
 * Show all lclasses as list, select if to be included or not, they can be maximized or minimized (or edit button).
 * 
 *
 */
public class MassListCreatorPanel extends JPanel 
{
	private static final long serialVersionUID = 1L;
	
	public MassListCreatorPanel()
	{
		JButton startAIButton = new JButton("Start");
  	startAIButton.setActionCommand("Start AI");
  	startAIButton.addActionListener(new ActionListener() {
		  public void actionPerformed(ActionEvent e) 
		  {
		  	startButton_actionPerformed(e);
		  }
	  });
  	this.add(startAIButton);
	}
	
	private void startButton_actionPerformed(ActionEvent e)
	{
		try 
		{
			ArrayList<AdductVO> allDefinedAdducts = (new AdductParser()).parse();
			ArrayList<LipidClassVO> allDefinedLipidClasses = (new LipidClassParser(allDefinedAdducts)).parse();
			MassListExporter exporter = new MassListExporter("D:\\Development\\eclipse_workspace\\LDA2_private\\massListCreation\\test2.xlsx", allDefinedLipidClasses);
			exporter.export();
			
		}
		catch (IOException | ChemicalFormulaException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
	}
	
}
