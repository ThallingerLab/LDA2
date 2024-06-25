package at.tugraz.genome.lda.masslist;

import java.io.FileOutputStream;
import java.io.IOException;

import javax.swing.JFrame;

import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.vos.AdductVO;

public class AdductExporter
{
	private AdductVO toExport_;
	
	public AdductExporter(AdductVO toExport)
	{
		this.toExport_ = toExport;
	}
	
	public void export()
	{
		try (FileOutputStream out= new FileOutputStream(buildAdductPath(toExport_.getFileName()));)
		{
			out.write("## The name of the adduct.\n".getBytes());
			out.write(String.format("%s=%s\n", AdductParser.PROPERTY_NAME, toExport_.getAdductName()).getBytes());
			out.write("## The chemical formula of the adduct. Chemical formulas prefixed with a '-' and '+' are subtracted from and added to the total mass, respectively. Syntax e.g. +Na-H2\n".getBytes());
			out.write(String.format("%s=%s\n", AdductParser.PROPERTY_FORMULA, toExport_.getFormulaString()).getBytes());
			out.write("## The charge state of the molecule with this adduct.\n".getBytes());
			out.write(String.format("%s=%s\n", AdductParser.PROPERTY_CHARGE, toExport_.getCharge()).getBytes());
		}
		catch (IOException ex) 
		{
			new WarningMessage(new JFrame(), "Error", "The export of the adduct definition file failed. Error message: "+ex.getMessage());
		}
	}
	
	public static String buildAdductPath(String fileName)
	{
		return AdductParser.ADDUCT_FOLDER+"/"+fileName;
	}
	
}
