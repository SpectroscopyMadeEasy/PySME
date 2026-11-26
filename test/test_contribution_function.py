"""Public spherical workflow integration test (performs stellar synthesis)."""
import numpy as np
from pysme.contribution_function import contribution_function
from pysme.sme_synth import SME_DLL


def check_spherical_workflow(sme, monkeypatch):
    """Shared assertions for the CI fixture and external frozen-model regression."""
    method_before = sme.atmo.method
    called=[]
    original=SME_DLL.GetLineDepressionContributionfunction
    def capture(self,*args,**kwargs):
        result=original(self,*args,**kwargs); called.append(result);return result
    monkeypatch.setattr(SME_DLL,'GetLineDepressionContributionfunction',capture)
    result=contribution_function(sme,kind='line_depression')[0]
    assert called
    np.testing.assert_allclose(result['contributions'],result['L']+result['Q']+result['N'])
    meta=result['native_metadata']
    np.testing.assert_allclose(meta[...,0],result['continuum_intensity'],rtol=1e-6,atol=1e-8)
    np.testing.assert_allclose(meta[...,1],result['intensity'],rtol=1e-6,atol=1e-8)
    np.testing.assert_allclose(result['contributions'].sum(-1)+result['boundary'],meta[...,0]-meta[...,1],rtol=1e-6,atol=1e-8)
    assert result['contributions'].shape[-1]==len(result['depth_edges'])-1
    assert sme.atmo.method == method_before  # Preparation used a private copy.
    emission=contribution_function(sme,kind='intensity')[0]
    np.testing.assert_allclose(emission['contributions'].sum(-1)+emission['boundary'],emission['native_metadata'][...,1],rtol=1e-10)
    return result, emission


def test_spherical_public_workflow(sme_2segments, monkeypatch):
    sme=sme_2segments
    sme.teff=4000; sme.logg=1; sme.vmic=2
    sme.atmo.source='marcs2012s_t2.0.sav'; sme.atmo.geom='SPH'
    sme.wave=[np.linspace(6554.,6555.,5)]
    sme.wran=[[6554.,6555.]]; sme.mu=[1.,.7]
    check_spherical_workflow(sme, monkeypatch)
