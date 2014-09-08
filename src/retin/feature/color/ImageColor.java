package retin.feature.color;

import retin.toolbox.core.*;
import retin.toolbox.document.*;

public class ImageColor
{
	static
	{
		Library.load("feature_color");
	}

	protected Parameters params;

	public ImageColor (String params) throws Exception
	{
		this.params = new Parameters(params);
	}

	public ImageColor (Parameters params)
	{
		this.params = params;
	}

	public native BytesMatrix runBytesMatrix (BytesMatrix input);
	public native FeatureMatrix runFeatureMatrix (FeatureMatrix input);

}
