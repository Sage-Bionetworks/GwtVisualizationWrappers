package org.gwtvisualizationwrappers.client.cytoscape;

/*
 * #%L
 * GwtCytoscapeJs
 * %%
 * Copyright (C) 2015 GwtCytoscapeJs
 * %%
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * #L%
 */

import com.google.gwt.dom.client.Element;
import com.google.gwt.user.client.ui.Widget;

public class CytoscapeGraph {
	/**
	 * Construct and show a cytoscape graph in the given container
	 * @param container
	 * @param params
	 */
	public CytoscapeGraph(Widget container, CytoscapeInitParams params) {
		show(container.getElement(), params);
	}
	private native void show(Element e, CytoscapeInitParams params) /*-{
		$wnd.jQuery(e).cytoscape(params);
	}-*/;
}
