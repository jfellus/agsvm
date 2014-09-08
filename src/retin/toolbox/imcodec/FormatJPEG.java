/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package retin.toolbox.imcodec;

import java.io.File;

/**
 *
 * @author romanegr
 */
public class FormatJPEG {
    public static void format (String file) throws Exception
	{
            if(file.toLowerCase().endsWith(".jpg"))
            {
                String filePpm = file.substring(0, file.length()-4) + ".ppm";
                String[] convJPEG2PPM = { "convert", file, filePpm};
                String[] convPPM2JPEG = { "convert", filePpm, file};
                String[] convRMJPEG = { "rm" , "-f", file };
                String[] convRMPPM = { "rm" , "-f", filePpm };

                if (new File(file).exists()) {
                    Runtime.getRuntime().exec(convJPEG2PPM);
                    if (new File(filePpm).exists()) {
                        Runtime.getRuntime().exec(convRMJPEG);
                        if(!new File(file).exists()) {
                            Runtime.getRuntime().exec(convPPM2JPEG);
                            if(new File(file).exists()) {
                                Runtime.getRuntime().exec(convRMPPM);
                                if(new File(filePpm).exists())
                                    throw new Exception ("Impossible de supprimer le fichier " + filePpm);
                            }
                            else
                                throw new Exception ("Impossible de créer le fichier " + file);
                        }
                        else
                            throw new Exception ("Impossible de supprimer le fichier " + file);
                    }
                    else
                        throw new Exception ("Impossible de créer le fichier " + filePpm);
                }
                else
                    throw new Exception ("Le fichier " + file + "n'existe pas");
            }
	}
}
