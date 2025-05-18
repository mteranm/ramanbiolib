import pytest
import numpy as np
from ramanbiolib.search import PeakMatchingSearch, SpectraSimilaritySearch

wavenumbers = np.array(list(range(450, 1800+1)))
pm_search = PeakMatchingSearch(wavenumbers)
spectra_search = SpectraSimilaritySearch(wavenumbers)

@pytest.mark.parametrize(
    "component",
    spectra_search.database['component'].unique()
)
def test_exact_component_spectra_search(component):
    spectrum = np.array(spectra_search.database[spectra_search.database['component'] == component]['intensity'].iloc[0])
    results = spectra_search.search(
        spectrum,
        unique_components_in_results=True,
        similarity_method="cosine_similarity",
    ).get_results(limit=None)
    
    # In case of ties consider first the component we know we are searching
    results['is_component'] = results['component'].apply(lambda x: 1 if x == component else 0)
    result = results.sort_values(["similarity_score", 'is_component'], ascending=False).iloc[0]

    assert result['component'] == component
    assert float("{:.2f}".format(result['similarity_score'])) == 1

@pytest.mark.parametrize(
    "component",
    pm_search.database['component'].unique()
)
def test_exact_component_pm_search(component):
    peaks = pm_search.database[pm_search.database['component'] == component]['peaks'].to_numpy()[0]
    peaks = list(filter(lambda x: x>=wavenumbers.min() and x<=wavenumbers.max(), peaks))
    results = pm_search.search(
        peaks,
        tolerance=5,
        class_filter=None,
        sort_score='PIUR',
        min_peak_intensity=0,
        tol_penalty="linear",
        unique_components_in_results=True
    ).get_results(limit=None)
    
    # In case of ties consider first the component we know we are searching
    results['is_component'] = results['component'].apply(lambda x: 1 if x == component else 0)
    result = results.sort_values(["PIUR", 'is_component'], ascending=False).iloc[0]

    assert result['component'] == component
    assert result['MR'] == 1
    assert result['RMR'] == 1
    assert result['IUR'] == 1