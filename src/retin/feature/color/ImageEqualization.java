package retin.feature.color;

import retin.toolbox.core.*;
import retin.toolbox.document.*;

public class ImageEqualization
{
	static
	{
		Library.load("feature_color");
	}

	protected Parameters params;

	public ImageEqualization (String params) throws Exception
	{
		this.params = new Parameters(params);
	}

	public ImageEqualization (Parameters params)
	{
		this.params = params;
	}

	public native BytesMatrix runBytesMatrix (BytesMatrix input);
	public native FeatureMatrix runFeatureMatrix (FeatureMatrix input);

}
