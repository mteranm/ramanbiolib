import time
import itertools
import pkg_resources

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from numpy.linalg import norm
from scipy.spatial import distance
from scipy import interpolate

class SearchResults:
    """
    Parent class for identification search results.

    Attributes
    ----------
    database : pandas.DataFrame
        The Raman spectra database dataframe
    wavenumbers : numpy.ndarray
        The numpy array with the specturm x-axis wavenumbers
    results_df : pandas.DataFrame
        The dataframe with the search results
    """
    def __init__(self, results_df, wavenumbers):
        """
        Parameters
        ----------
        results_df : pandas.DataFrame
            The dataframe with the search results
        wavenumbers : numpy.ndarray
            The numpy array with the specturm x-axis wavenumbers
        """
        self.database = pd.read_csv(
            pkg_resources.resource_filename(__name__, "db/raman_spectra_db.csv"),
            converters={
                "wavenumbers": lambda x: [float(v) for v in x.strip("[]").split(", ")],
                "intensity": lambda x: [float(v) for v in x.strip("[]").split(", ")]
            },
        )

        # Interpolate (if necessary)
        self.database['intensity'] = (
                self.database[['wavenumbers', 'intensity']]
                .apply(lambda x: list(interpolate.interp1d(x['wavenumbers'], x['intensity'])(wavenumbers)),
                    axis=1
                )
            )
        self.database['wavenumbers'] = self.database['wavenumbers'].apply(lambda x: list(wavenumbers))
        self.wavenumbers = wavenumbers
        self.results_df = results_df

    def _get_most_similar_type(self, k=5, sort_col='similarity_score'):
        """Gets the most similar type for the unkwnow spectrum.
        It uses k-NN with majority voting.

        Parameters
        ----------
        k : int
            The k nearest neighbours considered.
        sort_col : int
            The metric used to find the nearest neighbours.
        """
        top_k_results = self.results_df.drop_duplicates(['component'])[:k]
        top_k_results['type_lvl1'] = top_k_results['type'].apply(lambda x: x.split("/")[0])
        most_similar_classes = (
            top_k_results.groupby("type_lvl1").count()[['component']]
            .rename({"component": "count"}, axis=1)
            .sort_values('count', ascending=False)
            .reset_index()
            .merge(
                top_k_results
                .groupby("type_lvl1").mean()[[sort_col]]
                .sort_values(sort_col, ascending=False)
                .reset_index(),
                on=['type_lvl1']
            )
            .sort_values(['count', sort_col], ascending=False)
        )
        return most_similar_classes.iloc[0]
    
    def plot_results(self, y=[], n=3, peaks=None, height=600, sort_col='similarity_score'):
        """Plots the query spectrum/peaks and the top N results for comparison.

        Parameters
        ----------
        y : list
            The query spectrum intensity array
        n : int
            The N value to use the top N results
        peaks : list
            The query spectrum peaks position (in wavenumbers) list to plot
        height : int
            The plotly output plot height in pixels

        Raise
        ------
        Exception
            If neither y, nor peaks are empty/None.
        """
        if (y is None or len(y) == 0) and peaks is None:
            raise Exception("Either y or peaks must be NOT None")

        plot_df = pd.DataFrame({'wavenumbers': self.wavenumbers})
        if len(y) > 0:
            plot_df['search-query'] = y

        for i, row in self.results_df.sort_values(sort_col, ascending=False).reset_index(drop=True)[:n].iterrows():
            plot_df[f"{i+1}-{row['component']}"] = (
                self.database[(self.database['id'] == row['id'])]
                .iloc[0]
            )['intensity']
            
        fig = px.line(plot_df, x='wavenumbers', y=plot_df.drop("wavenumbers", axis=1).columns)
        if not peaks is None:
            for peak in peaks:
                fig = fig.add_vline(
                    x=peak, line_width=1, line_dash="dash", 
                    line_color="grey", 
                )
            fig = fig.add_trace(
                go.Scatter(
                    x=peaks,
                    y=(
                        y[np.searchsorted(self.wavenumbers, peaks)] if len(y) > 0 
                        else [1 for v in peaks]
                    ),
                    mode='markers',
                    name='peaks-query',
                    marker=dict(symbol='arrow-down', color='black'),
                )
            )
        fig.update_layout(
            xaxis_title="Wavenumbers (cm⁻¹)",
            yaxis_title="Intensity",
            legend_title="Component",
            height=height
        )
        return fig

    def plot_scores(self, height=500, metric='similarity_score'):
        """Plots the result scores by sorted index.

        Parameters
        ----------
        height : int
            The plotly output plot height in pixels
        metric : str
            The similarity score column name
        """
        return px.line(
            self.results_df, 
            y=metric, 
            hover_data=['id', 'component', 'laser', 'type', 'reference', 'source'],
            height=height
        )

class SpectraSearchResults(SearchResults):
    """
    Results object for spectra similarity search results.

    Attributes
    ----------
    database : pandas.DataFrame
        The Raman spectra database dataframe
    wavenumbers : numpy.ndarray
        The numpy array with the specturm x-axis wavenumbers
    results_df : pandas.DataFrame
        The dataframe with the search results
    query : pnumpy.ndarray
        The query spectrum intensity array
    """
    def __init__(self, results_df, query, wavenumbers):
        """
        Parameters
        ----------
        results_df : pandas.DataFrame
            The dataframe with the search results
        query : numpy.ndarray
            The query spectrum intensity array
        wavenumbers : numpy.ndarray
            The numpy array with the specturm x-axis wavenumbers
        """
        super().__init__(results_df, wavenumbers)
        self.query = query

    def get_results(self, limit=None, show_intensities=False):
        """Gets the search results dataframe.

        Parameters
        ----------
        limit : int
            The amount of results to show.
        show_intensities : bool
            The flag to indicate if the intensities column is returned, or not.
        """
        limit = len(self.results_df) if limit == None else limit
        results = self.results_df[:limit][['type', 'id', 'component', 'laser', 'reference', 'source', 'similarity_score', 'intensity']]
        return (
            results if show_intensities
            else results.drop("intensity", axis=1)
        )

    def get_most_similar_class(self, k=5):
        """Gets the most similar type for the unkwnow spectrum.
        It uses k-NN with majority voting.

        Parameters
        ----------
        k : int
            The k nearest neighbours considered.
        """
        return self._get_most_similar_type(k=k, sort_col='similarity_score')

    def plot_results(self, n=3, height=600, sort_col='similarity_score'):
        """Plots the query spectrum and the top N results for comparison.

        Parameters
        ----------
        n : int
            The N value to use the top N results
        height : int
            The plotly output plot height in pixels

        Raises
        ------
        Exception
            If query spectrum is empty/None.
        """
        return super().plot_results(y=self.query, n=n, height=height, sort_col=sort_col)

    def plot_scores(self, height=500):
        """Plots the result scores by sorted index.

        Parameters
        ----------
        height : int
            The plotly output plot height in pixels
        """
        return super().plot_scores(height=height)


class PMSearchResults(SearchResults):
    """
    Results object for spectra peak matching search results.

    Attributes
    ----------
    wavenumbers : numpy.ndarray
        The numpy array with the specturm x-axis wavenumbers
    results_df : pandas.DataFrame
        The dataframe with the search results
    peaks_a : numpy.ndarray
        The query peak positions (in wavenumbers) array
    assignments_df: pandas.DataFrame
        The dataframe containing the search result matching assignment for each db component
    """
    def __init__(self, results_df, wavenumbers, peaks_a, assignments_df):
        """
        Parameters
        ----------
        results_df : pandas.DataFrame
            The dataframe with the search results
        wavenumbers : numpy.ndarray
            The numpy array with the specturm x-axis wavenumbers
        peaks_a : numpy.ndarray
            The query peak positions (in wavenumbers) array
        assignments_df: pandas.DataFrame
            The dataframe containing the search result matching assignment for each db component
        """
        self.peaks_a = peaks_a
        self.assignments_df = assignments_df
        super().__init__(results_df, wavenumbers)

    def get_results(self, sort_col="IUR", limit=None):
        """Gets the search results dataframe.

        Parameters
        ----------
        sort_col : str
            The column used to rank the results [MR, RMR, IUR, PIUR]. Default IUR.
        limit : int
            The amount of results to show.
        """
        limit = len(self.results_df) if limit == None else limit
        results = self.results_df.copy()
        results = (
            results.sort_values(
                by=sort_col, ascending=False)[:limit]
                [[
                    'type', 'id', 'component', 'laser', 
                    'reference', 'source', 'MR', 'RMR', 'IUR', 'PIUR'
                ]]
        )
        return results

    def get_most_similar_class(self, k=5, sort_col='IUR'):
        """Gets the most similar type for the unkwnow spectrum.
        It uses k-NN with majority voting.

        Parameters
        ----------
        k : int
            The k nearest neighbours considered.
        sort_col : int
            The metric used to find the nearest neighbours.
        """
        return self._get_most_similar_type(k=k, sort_col=sort_col)

    def plot_pm_scores_matrix(
        self, limit=5, height=400, width=1000, color_continuous_scale="blues"
    ):
        """Plots the peaks matching scores matrix for the top N (limit) results.
        The matching scores are in range [0, 1] based on the matching distance and the penalization function defined in the search.

        Parameters
        ----------
        limit : int
            The amount of components to consider in the plot.
        height : int
            The plotly output plot height in pixels
        width : int
            The plotly output plot width in pixels
        color_continuous_scale : str
            The plotly output plot color_continuous_scale
        """
        wavenumbers_cols = list(
            filter(lambda col: str(col).isnumeric(), self.results_df.columns)
        )
        self.results_df['type_component'] = self.results_df.apply(
            lambda x: f"[{x['type']}] {x['component']}",
            axis=1
        ).to_numpy()
        return px.imshow(
            self.results_df[wavenumbers_cols].to_numpy()[:limit],
            x=[str(v) for v in wavenumbers_cols],
            y=self.results_df['type_component'][:limit],
            height=height,
            width=width,
            color_continuous_scale=color_continuous_scale,
        )

    def plot_results(self, n=3, height=600, query_spectrum=[], sort_col='IUR'):
        """Plots the query spectrum peaks positions and the top N results for comparison.

        Parameters
        ----------
        n : int
            The N value to use the top N results
        height : int
            The plotly output plot height in pixels
        query_spectrum : list, Optional
            The list of query spectrum intensity trace to be added to the comparison plot. It must match with db sectra wavenumbers dimension.

        Raise
        ------
        Exception
            If query spectrum and peaks are empty/None.
        """
        return super().plot_results(
            n=n, height=height, 
            peaks=self.peaks_a, y=np.array(query_spectrum),
            sort_col=sort_col
        )

    def plot_scores(self, height=500, metric='IUR'):
        """Plots the result scores by sorted index.

        Parameters
        ----------
        height : int
            The plotly output plot height in pixels
        metric : int
            The metric used in the plot. Values: MR, RMR, IUR, PIUR. Default IUR.
        """
        return super().plot_scores(height=height, metric='IUR')

    def get_assignments(self, component='default'):
        """Gets the query - db component peak matching assignments detail dataframe. 
        For each query peak position, that was matched in the db, the peak position found in the db component.

        Parameters
        ----------
        component : str
            The component name. If default, the top result component is used.
        """        
        component = self.results_df['component'].iloc[0] if component == 'default' else component
        # Get the component results
        component_assignment = (            
            self.assignments_df[
                self.assignments_df['component'] == component
            ][['id', 'component', 'laser', 'reference', 'source', 'peak_a', 'peak_b', 'diff']]
        )
        component_results = (            
            self.results_df[
                self.results_df['component'] == component
            ][['id', 'component', 'laser', 'reference', 'source', 'type']]
        )

        # Merge the best results to obtain the assignment of the best case for each component
        component_assignment = (
            component_assignment.merge(
                component_results,
                on=['id', 'component', 'laser', 'reference', 'source']
            )
        )
            
        return (
            component_assignment
            .rename({
                "peak_a": "peak_query",
                "peak_b": "peak_db",
                "diff": "distance"
            }, axis=1)
            .sort_values(["component", "laser", "peak_query"]).reset_index(drop=True)
        )

class _MatchingResult:
    """
    Peak matching result wrapper object.
    """
    def __init__(self, 
        no_pen_matching_scores, 
        no_pen_b_matching_scores, matching_scores, 
        b_matching_scores, matched_wavenumbers, 
        b_matched_wavenumbers
        ):
        self.no_pen_matching_scores = no_pen_matching_scores
        self.no_pen_b_matching_scores = no_pen_b_matching_scores
        self.matching_scores = matching_scores
        self.b_matching_scores = b_matching_scores
        self.matched_wavenumbers = matched_wavenumbers
        self.b_matched_wavenumbers = b_matched_wavenumbers

class ComponentSearch:
    """
    Parent class for component identification search.

    Attributes
    ----------
    database : pandas.DataFrame
        The database dataframe
    """
    def __init__(self, database):
        """
        Parameters
        ----------
        database : pandas.DataFrame
            The database dataframe
        """
        self.metadata = pd.read_csv(pkg_resources.resource_filename(__name__, "db/metadata_db.csv"))
        self.database = (
            database.merge(
                self.metadata[['id', 'source','reference', 'type', 'laser_wavelength']]
                .rename({"laser_wavelength": "laser"}, axis=1),
                on=['id']
            )
        )

    def list_db_components(self):
        """Lists the unique component values in the database"""
        return self.database['component'].unique()

    def list_db_types(self):
        """Lists the unique type values in the database"""
        return self.database['type'].unique()

    def list_db_lasers(self):
        """Lists the unique laser wavelenght values in the database"""
        return self.database['laser'].unique()

    def list_db_references(self):
        """Lists the unique article reference values in the database"""
        return self.database['reference'].unique()

class PeakMatchingSearch(ComponentSearch):
    """
    Implements a peak matching search over the peaks position database of Raman spectra for biomolecules.

    Attributes
    ----------
    database : pandas.DataFrame
        The peak positions database dataframe
    wavenumbers: numpy.ndarray
        The array of wavenumbers (integer) considered in the database. Max range 450 - 1800 cm⁻¹. Default value None, uses the default wavenumbers (range: 450 - 1800 cm⁻¹, step: 1 cm⁻¹)

    """
    def __init__(self, wavenumbers=None):
        """
        Parameters
        ----------
        wavenumbers: numpy.ndarray
            The array of wavenumbers (integer) considered in the database. Max range 450 - 1800 cm⁻¹. Default value None, uses the default wavenumbers (range: 450 - 1800 cm⁻¹, step: 1 cm⁻¹) 
        """
        database = pd.read_csv(
            pkg_resources.resource_filename(__name__, "db/raman_peaks_db.csv"),
            converters={
                "peaks": lambda x: [float(v) for v in x.strip("[]").split(", ")],
                "intensity": lambda x: [float(v) for v in x.strip("[]").split(", ")],
            }
        )
        self.max_peaks_len = database['peaks'].str.len().max()
        self.wavenumbers = wavenumbers if not wavenumbers is None else np.arange(450, 1801)
        self.__MAX_ASSIGNMENT_COMBINATIONS = 5000
        super().__init__(database)

    def __get_intersection_union_score(
        self,
        component,
        peaks_a,
        peaks_b,
        matching_scores,
        b_matching_scores,
        matched_wavenumbers,
        b_matched_wavenumbers,
    ):
        union = np.union1d(peaks_a, peaks_b)
        intersection = np.union1d(
            matched_wavenumbers, b_matched_wavenumbers
        )  # The union of the matched wavenumber in both sides is the intersection considering a tolerance
        matching_scores_merged = {**matching_scores, **b_matching_scores}
        intersection_scores = [
            matching_scores_merged[wavenumber] for wavenumber in intersection
        ]
        return sum(intersection_scores) / len(union)

    def __get_matching_score(self, peak_a, peak_b, tolerance, tol_penalty=None):
        # No tolerance penalty
        if tol_penalty == None:
            return 1
        # Tolerance penalty calculation
        diff = np.abs(peak_a - peak_b)
        if tol_penalty == "linear":
            return np.clip(1 - (diff / (tolerance + 1)), a_min=0, a_max=1)
        elif tol_penalty == "inverse_power":
            return np.clip(
                ((tolerance + 1) - diff) / ((tolerance + 1) * (1 + diff)),
                a_min=0,
                a_max=1,
            )
        else:
            raise Exception(f"Unsupported tolerance method {tol_penalty}")

    def __get_assignation_results(self, part_df, increasing_assingation_a, assigned_peaks_b):
        # Get best assignation dataframe
        assign_merged = part_df.merge(
            pd.DataFrame({
                "peak_a": increasing_assingation_a, 
                "peak_b": assigned_peaks_b
            }),
            on=['peak_a', 'peak_b']
        )

        return {
            "a_matched": len(assign_merged['peak_a']),
            "score_sum": assign_merged['score_pen'].sum(),
            "peaks_a": assign_merged['peak_a'].tolist(),
            "peaks_b": assign_merged['peak_b'].tolist(),
            "scores": assign_merged['score'].tolist(),
            "scores_pen": assign_merged['score_pen'].tolist()
        }

    def __get_peak_a_assignation_combinations(self, matching_df_grouped):
        matching_df_grouped_reindex = matching_df_grouped.reset_index()
        peak_b_duplicates = matching_df_grouped_reindex[matching_df_grouped_reindex['len_peaks_a'] > 1]

        # Create 0 to max peak_b duplicate integer number combinations
        # only for peak_b duplicated cases 
        # (this is more efficient than creating combinations for duplicated and not duplicated)
        b_duplicated_assignation_combinations = list(
            itertools.product(
                list(range(peak_b_duplicates['len_peaks_a'].max())), 
                repeat=len(peak_b_duplicates)
            )
        )

        if len(b_duplicated_assignation_combinations) > self.__MAX_ASSIGNMENT_COMBINATIONS:
            b_duplicated_assignation_combinations[:self.__MAX_ASSIGNMENT_COMBINATIONS]
            print(f"Warning! Too many duplicated assignments combinations. Not considering all cases for deduplication. Reduce the tolerance or the amount of close peaks.")

        # As combinations are created considering the max duplicates for a single peak, we need to filter
        # the cases considering only the maxium values for each peak duplicates amount
        groups_len = peak_b_duplicates['len_peaks_a'].tolist()
        b_assig_comb_filtered = list(
            filter(lambda x: not False in [ item < groups_len[i] for i, item in enumerate(x)], 
            b_duplicated_assignation_combinations)
        )

        # Get full combinations considering peak_b duplicated cases and not duplicated cases
        default_assignation_comb = [
            matching_df_grouped_reindex.iloc[index]['peak_a'][0]
            for index in matching_df_grouped_reindex.index.tolist()
        ]
        assignation_combs  = []
        for comb in b_assig_comb_filtered:
            new_full_combination = default_assignation_comb.copy()
            for i, position in enumerate(comb):
                new_full_combination[peak_b_duplicates.index[i]] = peak_b_duplicates.iloc[i]['peak_a'][position]
            assignation_combs.append(new_full_combination)

        # Get the assignations that are maximizing the len of unique assigned peak_a
        max_matching = np.max([len(np.unique(comb)) for comb in assignation_combs])
        most_matching_assingations = list(filter(lambda c: len(np.unique(c)) == max_matching, assignation_combs))

        # Fiter to only assignations with increasing wavenumbers
        increasing_assingations_a = list(filter(lambda c: np.all(np.diff(c) >= 0), most_matching_assingations))

        return increasing_assingations_a

    def __part_assignation_deduplication(self, part_df):
        # Group peak_a assignations by peak b
        matching_df_grouped = pd.DataFrame(part_df.groupby("peak_b")['peak_a'].apply(list))
        matching_df_grouped['len_peaks_a'] = matching_df_grouped['peak_a'].apply(lambda x: len(x))

        # Create peak_a assignations combinations
        increasing_assingations_a = self.__get_peak_a_assignation_combinations(matching_df_grouped)

        # Calculate assignation score for each combination
        assignations_df = pd.DataFrame([
            self.__get_assignation_results(part_df, increasing_assingation_a, matching_df_grouped.index)
            for increasing_assingation_a in increasing_assingations_a
        ])

        # Get best assignation
        return assignations_df.sort_values(by=['a_matched', 'score_sum'], ascending=False).iloc[0]

    def __peak_assignations_deduplication(self, matching_df, tolerance):
        # Group peaks a by close groups
        grouped_peak_a = []
        previous = 0
        current_group = []
        assignments = []
        for peak_a in matching_df['peak_a'].unique():
            if peak_a - previous <= tolerance*2:
                current_group.append(previous)
            elif len(current_group) > 0:
                current_group.append(previous)
                grouped_peak_a.append(current_group)
                current_group = []
            previous = peak_a

        ## Split df
        parts = []
        last_border_peak_a = 0
        for i, group in enumerate(grouped_peak_a):
            if i == 0:
                continue
            parts.append(matching_df[(matching_df['peak_a'] >= last_border_peak_a)&(matching_df['peak_a'] < group[0])])
            last_border_peak_a = group[0]

        last_part = matching_df[matching_df['peak_a'] >= last_border_peak_a]
        if len(last_part) > 0:
            parts.append(last_part)

        matched_wavenumbers = []
        b_matched_wavenumbers = []
        a_matching_scores = {}
        b_matching_scores = {}
        a_matching_scores_pen = {}
        b_matching_scores_pen = {}
        for part_df in parts:
            peak_b_duplicates = (len(part_df['peak_b'].unique()) < len(part_df['peak_b']))
            if peak_b_duplicates:
                best_assignation = self.__part_assignation_deduplication(part_df)
                assignments.append(best_assignation)
                matched_wavenumbers.extend([int(v) for v in best_assignation['peaks_a']])
                b_matched_wavenumbers.extend([int(v) for v in best_assignation['peaks_b']])
                for i, score in enumerate(best_assignation['scores']):
                    a_matching_scores[int(best_assignation['peaks_a'][i])] = score
                for i, score in enumerate(best_assignation['scores']):
                    b_matching_scores[int(best_assignation['peaks_b'][i])] = score
                for i, score in enumerate(best_assignation['scores_pen']):
                    a_matching_scores_pen[int(best_assignation['peaks_a'][i])] = score
                for i, score in enumerate(best_assignation['scores_pen']):
                    b_matching_scores_pen[int(best_assignation['peaks_b'][i])] = score
            else:
                assignments.append(part_df)
                for peak_a, score in  part_df[['peak_a', 'score']].to_numpy():
                    a_matching_scores[peak_a] = score
                for peak_b, score in  part_df[['peak_b', 'score']].to_numpy():
                    b_matching_scores[peak_b] = score
                for peak_a, score in  part_df[['peak_a', 'score_pen']].to_numpy():
                    a_matching_scores_pen[peak_a] = score
                for peak_b, score in  part_df[['peak_b', 'score_pen']].to_numpy():
                    b_matching_scores_pen[peak_b] = score
                matched_wavenumbers.extend(part_df['peak_a'].tolist())
                b_matched_wavenumbers.extend(part_df['peak_b'].tolist())

        return (
            a_matching_scores,
            b_matching_scores,
            a_matching_scores_pen,
            b_matching_scores_pen,
            matched_wavenumbers,
            b_matched_wavenumbers,
            pd.concat(assignments)
        )
     
    def __get_component_matched_peaks(
        self, peaks_a, peaks_b, tolerance, tol_penalty="linear"
    ):
        matching_pairs =  []
        for i, peak_a in enumerate(peaks_a):
            # Tolerance peak matching search
            for peak_b in peaks_b:
                if peak_b <= peak_a + tolerance and peak_b >= peak_a - tolerance:
                    diff = np.abs(peak_a - peak_b)
                    matching_score = self.__get_matching_score(
                        peak_a, peak_b, tolerance, tol_penalty=None
                    )
                    matching_score_pen = self.__get_matching_score(
                        peak_a, peak_b, tolerance, tol_penalty=tol_penalty
                    )
                    matching_pairs.append({
                        "peak_a": peak_a, "peak_b": peak_b, "diff": diff, 
                        "score": matching_score, "score_pen": matching_score_pen
                    })

        # Return empty if there is no matching
        if len(matching_pairs) == 0:
            return _MatchingResult({}, {}, {}, {}, [], []), pd.DataFrame()

        # Fix repeated matching (single peak_b matched to multiple peak_a)
        matching_pairs_df = pd.DataFrame(matching_pairs)
        peak_b_duplicates = (len(matching_pairs_df['peak_b'].unique()) < len(matching_pairs_df['peak_b']))
        if tolerance > 0 and peak_b_duplicates:
            a_matching_scores, b_matching_scores, a_matching_scores_pen, b_matching_scores_pen, matched_wavenumbers, b_matched_wavenumbers, assignments_df = (
                self.__peak_assignations_deduplication(matching_pairs_df, tolerance)
            )
        else:
            a_matching_scores = { peak_a:score for peak_a, score in matching_pairs_df[['peak_a', 'score']].to_numpy() }
            b_matching_scores = { peak_b:score for peak_b, score in matching_pairs_df[['peak_b', 'score']].to_numpy() }
            a_matching_scores_pen = { peak_a:score for peak_a, score in matching_pairs_df[['peak_a', 'score_pen']].to_numpy() }
            b_matching_scores_pen = { peak_b:score for peak_b, score in matching_pairs_df[['peak_b', 'score_pen']].to_numpy() }
            matched_wavenumbers = matching_pairs_df['peak_a'].tolist()
            b_matched_wavenumbers = matching_pairs_df['peak_b'].tolist()
            assignments_df = matching_pairs_df

        return _MatchingResult(
            { peak_a: a_matching_scores[peak_a] if peak_a in a_matching_scores else 0 for peak_a in peaks_a },
            { peak_b: b_matching_scores[peak_b] if peak_b in b_matching_scores else 0 for peak_b in peaks_b },
            { peak_a: a_matching_scores_pen[peak_a] if peak_a in a_matching_scores_pen else 0 for peak_a in peaks_a },
            { peak_b: b_matching_scores_pen[peak_b] if peak_b in b_matching_scores_pen else 0 for peak_b in peaks_b },
            matched_wavenumbers,
            b_matched_wavenumbers
        ), assignments_df

    def search(
        self,
        peaks_a,
        tolerance=0,
        class_filter=None,
        sort_score='IUR',
        min_peak_intensity=0,
        tol_penalty="linear",
        unique_components_in_results=True
    ):
        """Finds the most similar component in the database using a peak matching score.

        Considering query spectrum peaks Pa and DB spectrum peaks Pb.
        Metrics:

        - MR (Matching Ratio) = intersection(Pa, Pb)/len(Pa)
        - RMR (Reverse Matching Ratio) = intersection(Pa, Pb)/len(Pb)
        - IUR (Intersection Union Ratio) = intersection(Pa, Pb)/union(Pa, Pb)
        - PIUR (Penalized Intersection Union Ratio) = penalized_intersection(Pa, Pb)/union(Pa, Pb)

        Parameters
        ----------
        peaks_a : list
            The list of unknown spectrum peaks positions in wavenumbers (integer)
        tolerance: int
            The simmetrical maximum distance tolerance for peak matching. Default 0, only exact match.
        class_filter: list
            The list of classes to consider in the results. Default value is None, all cases considered.
        sort_score: str
            The matching score metric to use for the results ranking. Values: MR, RMR, IUR, PIUR. Default IUR.
        min_peak_intensity: float
            The minimum intensity value for database peaks considering the the matching (between 0-1). Default value 0.
        tol_penalty: str
            The type of penalty function applied for PIUR calculation. Values: 'linear' or 'inverse_power'. Default linear.
        unique_components_in_results: bool
            The flag to indicate if the results have duplicates measures results for each component. True when only the best result for each component, in the database, is shown in the results. False otherwise. Default True.
        """
        # Peaks size check
        if len(peaks_a) > 2*self.max_peaks_len:
            print(f"Warning! Too many peaks [{len(peaks_a)}]. The maximum peaks length in the db is {self.max_peaks_len}. Consider to reduce the number of peaks in the search query.")

        if len(peaks_a) > 3*self.max_peaks_len:
            raise Exception(f"Too many peaks [{len(peaks_a)}]. The maximum peaks length in the db is {self.max_peaks_len}. Reduce the number of peaks in the search query.")

        # Filter components by class (if configured)
        filtered_db = (
            self.database[self.database["type"].str.contains('|'.join(class_filter))]
            if class_filter != None
            else self.database
        )

        # Components search
        peaks_similarities = []
        assignments = []

        for i, row in filtered_db.iterrows():
            component = row['component']
            peaks_b = row["peaks"]
            peaks_b = list(
                filter(
                    lambda x: (x >= self.wavenumbers.min())
                    & (x <= self.wavenumbers.max()),
                    peaks_b,
                )
            )
            if len(peaks_b) == 0:
                continue

            matching_results, c_assignments_df = self.__get_component_matched_peaks(
                peaks_a, peaks_b, tolerance, tol_penalty=tol_penalty
            )

            row_dict = { int(peak): score for peak, score in matching_results.matching_scores.items() }
            row_dict['id'] = float(row['id'])
            row_dict["component"] = component
            row_dict["laser"] = row['laser']
            row_dict["reference"] = row['reference']
            row_dict["source"] = row['source']
            row_dict["type"] = row["type"]
            row_dict["matched"] = len(np.unique(matching_results.matched_wavenumbers))
            row_dict["MR"] = row_dict["matched"] / len(peaks_a)
            row_dict["RMR"] = row_dict["matched"] / len(peaks_b)
            row_dict["a_not_matched"] = len(
                np.unique(list(filter(lambda x: x not in matching_results.matched_wavenumbers, peaks_a)))
            )
            row_dict["b_not_matched"] = len(
                list(filter(lambda x: x not in matching_results.b_matched_wavenumbers, peaks_b))
            )
            row_dict["peaks_b"] = peaks_b
            row_dict["peaks_a"] = peaks_a
            row_dict["b_matched_wavenumbers"] = matching_results.b_matched_wavenumbers
            row_dict["a_matched_wavenumbers"] = matching_results.matched_wavenumbers
            row_dict["a_matching_scores"] = matching_results.matching_scores
            row_dict["b_matching_scores"] = matching_results.b_matching_scores
            row_dict["PIUR"] = self.__get_intersection_union_score(
                component,
                peaks_a,
                peaks_b,
                matching_results.matching_scores,
                matching_results.b_matching_scores,
                matching_results.matched_wavenumbers,
                matching_results.b_matched_wavenumbers,
            )
            row_dict["IUR"] = self.__get_intersection_union_score(
                component,
                peaks_a,
                peaks_b,
                matching_results.no_pen_matching_scores,
                matching_results.no_pen_b_matching_scores,
                matching_results.matched_wavenumbers,
                matching_results.b_matched_wavenumbers,
            )
            peaks_similarities.append(row_dict)

            c_assignments_df['id'] = float(row['id'])
            c_assignments_df['component'] = component
            c_assignments_df["laser"] = row['laser']
            c_assignments_df["reference"] = row['reference']
            c_assignments_df["source"] = row['source']
            assignments.append(c_assignments_df)

        # Results formatting
        detailed_df = pd.DataFrame(
            peaks_similarities,
        )

        # Get the max similarity score for each component results
        # As we can have different measurements of the same components we can have
        # multiple similarity for a single component
        if unique_components_in_results:
            max_score_df = (
                detailed_df.groupby(["component"])
                .max(sort_score).reset_index()[['component', sort_score]]
            )
            detailed_df =(
                detailed_df.merge(max_score_df, on=['component', sort_score])
                # If there are multiple cases with the same max score just drop duplicates
                .drop_duplicates(subset=['component'])
            )

        detailed_df = detailed_df.sort_values(
            by=[sort_score, "matched"], ascending=False
        )

        return PMSearchResults(
            detailed_df.reset_index(drop=True), 
            self.wavenumbers,
            peaks_a,
            pd.concat(assignments)
        )

class SpectraSimilaritySearch(ComponentSearch):
    """
    Implements a spectra similarity search over the database of Raman spectra for biomolecules.

    Attributes
    ----------
    database : pandas.DataFrame
        The Raman spectra database dataframe
    wavenumbers: numpy.ndarray
        The array of wavenumbers values for the spectra traces. Max range 450 - 1800 cm⁻¹. If a custom value is set, the db spectra are linearly interpolated to the new wavenumbers list. Default value None, uses the default wavenumbers (range: 450 - 1800 cm⁻¹, step: 1 cm⁻¹) 
    """
    def __init__(self, wavenumbers=None):
        """
        Parameters
        ----------
        wavenumbers: numpy.ndarray, Optional
            The array of wavenumbers values for the spectra traces. Max range 450 - 1800 cm⁻¹. If a custom value is set, the db spectra are linearly interpolated to the new wavenumbers list. Default value None, uses the default wavenumbers (range: 450 - 1800 cm⁻¹, step: 1 cm⁻¹) 
        """
        database = pd.read_csv(
            pkg_resources.resource_filename(__name__, "db/raman_spectra_db.csv"),
            converters={
                "wavenumbers": lambda x: [float(v) for v in x.strip("[]").split(", ")],
                "intensity": lambda x: [float(v) for v in x.strip("[]").split(", ")]
            },
        )

        if not wavenumbers is None:
            database['intensity'] = (
                database[['wavenumbers', 'intensity']]
                .apply(lambda x: list(interpolate.interp1d(x['wavenumbers'], x['intensity'])(wavenumbers)),
                    axis=1
                )
            )
            database['wavenumbers'] = database['wavenumbers'].apply(lambda x: list(wavenumbers))
        
        self.wavenumbers = np.array(database.iloc[0]['wavenumbers'])
        super().__init__(database)

    def cosine_similarity(self, spectra_a, spectra_b):
        """Gets the cosine similarity between two spectra. 
        Parameters
        ----------
        spectra_a: numpy.ndarray
            The first spectrum to find the similarity
        spectra_b: numpy.ndarray
            The second spectrum to find the similarity

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Cosine_similarity
        """
        return np.dot(spectra_a,spectra_b)/(norm(spectra_a)*norm(spectra_b))
    
    def slk_similarity(self, x, z, w=25):
        """Gets the spectral linear kernel (SLK) similarity between two spectra.  

        Parameters
        ----------
        x: numpy.ndarray
            The first spectrum to find the similarity
        z: numpy.ndarray
            The second spectrum to find the similarity
        w: int
            The window parameter value (The final window size is 2*w+1 around each point). Default 25.

        References
        ----------
        ..  [2] Khan et al., 'New similarity metrics for Raman spectroscopy', Chemom. Intell. Lab. Syst., vol. 114, pp. 99–108, May 2012, doi: 10.1016/j.chemolab.2012.03.007.
        """        
        # We need to pad 0 on the edges to be able to have symmetrical 2*W windows for every elemenent.
        x_pad = np.pad(x, (w, w), 'constant', constant_values=(0, 0))
        z_pad = np.pad(z, (w, w), 'constant', constant_values=(0, 0))

        # Create the 2*W+1 windows for each point of the spectrum (i-W, i+W)
        x_ws = np.lib.stride_tricks.sliding_window_view(x_pad, w*2+1)
        z_ws = np.lib.stride_tricks.sliding_window_view(z_pad, w*2+1)

        # Calculate the similarity for each point
        # x[i]*z[i] + (sum((x[i]-x[j])*z[i]-z[j])) over the window)
        # Calculate (x[i]-x[j])*(z[i]-z[j]) vectorized is much faster
        point_similarity = np.sum(((x_ws[:, w] - x_ws.transpose())*(z_ws[:, w] - z_ws.transpose())), axis=0) + (x*z)
        
        # Calculate the final score
        return np.sum(point_similarity)/len(x)

    def __calculate_similarity(self, method, spectra_a, spectra_b, similarity_params):
        if method == 'euclidean':
            return distance.euclidean(spectra_a, spectra_b)
        elif method == 'rmse':
            return np.sqrt(np.mean(np.square(np.abs(spectra_a - spectra_b))))
        elif method == 'slk':
            spectra_a = np.array(spectra_a)
            spectra_b = np.array(spectra_b)
            return (
                self.slk_similarity(spectra_a, spectra_b, w=similarity_params)/
                np.sqrt(
                    self.slk_similarity(spectra_a, spectra_a, w=similarity_params)* 
                    self.slk_similarity(spectra_b, spectra_b, w=similarity_params)
                )
            )
        elif method == 'cosine_similarity':
            return self.cosine_similarity(spectra_a, spectra_b)

        raise Exception(f"Wrong similarity method value: {method}")

    def search(
        self,
        spectra_a,
        class_filter=None,
        unique_components_in_results=True,
        similarity_method="slk",
        similarity_params=50
        ):
        """Finds the most similar component in the database using a spectra similarity score.

        Parameters
        ----------
        spectra_a : list
            The unknown spectrum intensity trace. The values must be interpolated to match the wavenumbers list specified when creating the SpectraSimilaritySearch object.
        class_filter: list
            The list of classes to consider in the results. Default value is None, all cases considered.
        unique_components_in_results: bool
            The flag to indicate if the results have duplicates measures results for each component. True when only the best result for each component, in the database, is shown in the results. False otherwise. Default True.
        similarity_method: str
            The similarity score used for results ranking. Values: 'euclidean', 'cosine_similarity', 'slk'. Default: 'slk'.
        similarity_params: 
             In slk case, the value of the window (w) parameter. Is ignored in other cases.
       """     
        # Filter components by class (if configured)
        filtered_db = (
            self.database[self.database["type"].str.contains('|'.join(class_filter))]
            if class_filter != None
            else self.database
        )

        # Init results
        search_results = filtered_db[
            ['id', 'component', 'intensity', 'laser', 'reference', 'source', 'type']
        ].copy()

        # Min-max normalization over spectra_a to fit with db normalization
        spectra_a = (spectra_a - np.min(spectra_a))/(np.max(spectra_a) - np.min(spectra_a))

        # Calculate the similarity score for each database item
        search_results['similarity_score'] = (
            search_results['intensity']
            .apply(lambda spectra_b: self.__calculate_similarity(similarity_method, spectra_a, spectra_b, similarity_params))
        )

        # Get the max similarity score for each component results
        # As we can have different measurements of the same components we can have
        # multiple similarity for a single component
        if unique_components_in_results:
            max_score_df = (
                search_results.groupby(["component"])
                .max("similarity_score").reset_index()[['component', 'similarity_score']]
            )
            search_results = (
                search_results.merge(max_score_df, on=['component', 'similarity_score'])
                # If there are multiple cases with the same max score just drop duplicates
                .drop_duplicates(subset=['component'])
            )

        return SpectraSearchResults(
            search_results.sort_values(by='similarity_score', ascending=False).reset_index(drop=True),
            spectra_a,
            self.wavenumbers
        )