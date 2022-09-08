# -*- coding: utf-8 -*-

import plotly.graph_objects as go


import spectral_match

def spectra_comparison(
    sample_mz: list,
    sample_sig: list,
    ref_mz: list,
    ref_sig: list,
    sample_label: str = 'Compound 1',
    ref_label: str = 'Compound 2',
    html_out: bool = False):
    """
    

    Parameters
    ----------
    
    Mass Spectrometry Common Outputs
    --------------------------------
    These are all common outputs of exporting data from mass spec. instruments
    MassHunter, Chromeleon, and NIST library software all export in the format
        
        sample_mz : list
            list of all detected m/z values detected in sample mass spectrum
        sample_sig : list
            list of all signals corresponding to m/z in sample mass spectrum
        ref_mz : list
            list of all detected m/z values detected in reference mass spectrum
        ref_sig : list
            list of all signals corresponding to m/z in reference mass spectrum
    
    sample_label : str, optional
        Label for sample compound
        The default is 'Compound 1'.
    ref_label : str, optional
        Label for reference compound 
        The default is 'Compound 2'.
    html_out : bool, optional
        Toggle for creating the plot as interactive html file.
        The default is False and will not create a file.

    Returns
    -------
    fig : plotly.graph_object
        graph object created in plotly

    """
    
    norm_sample_sig = vector_norm(sample_sig)
    norm_ref_sig = vector_norm(ref_sig)
    
    # properly orients reference spectrum downwards
    down_norm_ref_sig = [-1 * value for value in norm_ref_sig]
    
    #match factor calculations
    match_factor = round(match(sample_sig,sample_mz,ref_sig,ref_mz),3)
    fwd_match_factor = round(fwd_match(sample_sig,sample_mz,ref_sig,ref_mz),3)
    rev_match_factor = round(rev_match(sample_sig,sample_mz,ref_sig,ref_mz),3)
    
    # plotly visualization
    fig = go.Figure()
    fig.add_trace(go.Bar(name = sample_label, x = sample_mz, y = norm_sample_sig, offset = 0))
    fig.add_trace(go.Bar(name = ref_label, x = ref_mz, y = down_norm_ref_sig , offset = 0))
    
    fig.update_layout(hovermode="x")
    fig.update_layout(
        title = f"<b>{sample_label} vs. {ref_label} Mass Spectra Comparison</b>",
        title_x = 0.5,
        xaxis_title = '<br>'.join([
             '<b>m/z-value</b>',
             f'Match Factor: {match_factor}',
             f'Reverse Match Factor: {rev_match_factor}',
             f'Forward Match Fator: {fwd_match_factor}']),
        yaxis_title="<b>Normalized Intensity</b>",
        legend_title="<b>Legend Title<b>"
    )
    
    # save plot as interactive html file
    if html_out == True:
        fig.write_html(f"{sample_label}_vs_{ref_label}_ms_comparison.html")
    fig.show()
    return fig
               


