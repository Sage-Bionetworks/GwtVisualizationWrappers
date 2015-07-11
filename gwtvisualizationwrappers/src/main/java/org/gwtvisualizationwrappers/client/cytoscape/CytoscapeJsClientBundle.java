package org.gwtvisualizationwrappers.client.cytoscape;

/*
 * #%L
 * GwtVisualizationWrapper
 * %%
 * Copyright (C) 2015 - 2016 GwtVisualizationWrapper
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


import com.google.gwt.core.client.GWT;
import com.google.gwt.resources.client.ClientBundle;
import com.google.gwt.resources.client.TextResource;

/**
 * @author Jay Hodgson
 */
public interface CytoscapeJsClientBundle extends ClientBundle {

    static final CytoscapeJsClientBundle INSTANCE = GWT.create(CytoscapeJsClientBundle.class);

    @Source("resource/js/2.4.2-cytoscape.closurecompiled.js")
    TextResource cytoscape2_4_2();

    @Source("resource/js/jquery-1.11.2.min.cache.js")
    TextResource jQuery();
    
    @Source("resource/js/sage-cytoscape-utils.js")
    TextResource sageCytoscapeUtils();
}
