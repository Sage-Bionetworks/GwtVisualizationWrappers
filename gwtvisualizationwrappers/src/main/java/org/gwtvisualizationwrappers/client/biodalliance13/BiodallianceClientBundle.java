package org.gwtvisualizationwrappers.client.biodalliance13;

import com.google.gwt.core.client.GWT;
import com.google.gwt.resources.client.ClientBundle;
import com.google.gwt.resources.client.TextResource;

/**
 * @author Jay Hodgson
 */
public interface BiodallianceClientBundle extends ClientBundle {

    static final BiodallianceClientBundle INSTANCE = GWT.create(BiodallianceClientBundle.class);
    
    @Source("resource/js/0.13.dev-dalliance-compiled.js")
    TextResource biodalliance0_13dev();
}
