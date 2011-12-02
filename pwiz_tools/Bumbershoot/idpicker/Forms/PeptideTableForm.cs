﻿//
// $Id$
//
// The contents of this file are subject to the Mozilla Public License
// Version 1.1 (the "License"); you may not use this file except in
// compliance with the License. You may obtain a copy of the License at
// http://www.mozilla.org/MPL/
//
// Software distributed under the License is distributed on an "AS IS"
// basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
// License for the specific language governing rights and limitations
// under the License.
//
// The Original Code is the IDPicker project.
//
// The Initial Developer of the Original Code is Matt Chambers.
//
// Copyright 2010 Vanderbilt University
//
// Contributor(s): Surendra Dasari
//

using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.Threading;
using DigitalRune.Windows.Docking;
using NHibernate.Linq;
using BrightIdeasSoftware;
using PopupControl;
using IDPicker.DataModel;
using IDPicker.Controls;

namespace IDPicker.Forms
{
    public partial class PeptideTableForm : BaseTableForm
    {
        #region Wrapper classes for encapsulating query results

        public class AggregateRow : Row
        {
            public int Spectra { get; private set; }
            public int DistinctMatches { get; private set; }
            public int DistinctPeptides { get; private set; }
            public int Proteins { get; private set; }
            public string ProteinAccessions { get; private set; }
            public string ProteinGroups { get; private set; }

            public static int ColumnCount = 6;
            public static string Selection = "SELECT " +
                                             "COUNT(DISTINCT psm.Spectrum.id), " +
                                             "COUNT(DISTINCT psm.DistinctMatchKey), " +
                                             "COUNT(DISTINCT psm.Peptide.id), " +
                                             "COUNT(DISTINCT pro.id), " +
                                             "DISTINCT_GROUP_CONCAT(pro.Accession), " +
                                             "DISTINCT_GROUP_CONCAT(pro.ProteinGroup)";

            #region Constructor
            public AggregateRow (object[] queryRow, DataFilter dataFilter)
            {
                int column = -1;
                Spectra = Convert.ToInt32(queryRow[++column]);
                DistinctMatches = Convert.ToInt32(queryRow[++column]);
                DistinctPeptides = Convert.ToInt32(queryRow[++column]);
                Proteins = Convert.ToInt32(queryRow[++column]);
                ProteinAccessions = Convert.ToString(queryRow[++column]);
                ProteinGroups = Convert.ToString(queryRow[++column]);
                DataFilter = dataFilter;
            }
            #endregion
        }

        public class PeptideGroupRow : AggregateRow
        {
            public int PeptideGroup { get; private set; }

            #region Constructor
            public PeptideGroupRow (object[] queryRow, DataFilter dataFilter)
                : base(queryRow, dataFilter)
            {
                int column = AggregateRow.ColumnCount - 1;
                PeptideGroup = Convert.ToInt32(queryRow[++column]);
            }
            #endregion
        }

        public class DistinctPeptideRow : AggregateRow
        {
            public Peptide Peptide { get; private set; }

            #region Constructor
            public DistinctPeptideRow (object[] queryRow, DataFilter dataFilter)
                : base(queryRow, dataFilter)
            {
                int column = AggregateRow.ColumnCount - 1;
                Peptide = (DataModel.Peptide) queryRow[++column];
            }
            #endregion
        }

        public class DistinctMatchRow : AggregateRow
        {
            public Peptide Peptide { get; private set; }
            public DistinctMatchKey DistinctMatch { get; private set; }
            public PeptideSpectrumMatch PeptideSpectrumMatch { get; private set; }

            #region Constructor
            public DistinctMatchRow (object[] queryRow, DataFilter dataFilter)
                : base(queryRow, dataFilter)
            {
                int column = AggregateRow.ColumnCount - 1;
                Peptide = (Peptide) queryRow[++column];
                PeptideSpectrumMatch = (PeptideSpectrumMatch) queryRow[++column];
                DistinctMatch = new DistinctMatchKey(Peptide, PeptideSpectrumMatch,
                                                     dataFilter.DistinctMatchFormat,
                                                     (string) queryRow[++column]);
            }

            #endregion
        }

        public class PeptideInstanceRow : Row
        {
            public PeptideInstance PeptideInstance { get; private set; }

            #region Constructor
            public PeptideInstanceRow (object queryRow, DataFilter dataFilter)
            {
                PeptideInstance = (DataModel.PeptideInstance) queryRow;
            }
            #endregion
        }

        class PivotData
        {
            public int Spectra { get; private set; }
            public int DistinctMatches { get; private set; }
            public int DistinctPeptides { get; private set; }

            #region Constructor
            public PivotData () { }
            public PivotData (object[] queryRow)
            {
                Spectra = Convert.ToInt32(queryRow[3]);
                DistinctMatches = Convert.ToInt32(queryRow[4]);
                DistinctPeptides = Convert.ToInt32(queryRow[5]);
            }
            #endregion
        }

        struct TotalCounts
        {
            public int PeptideGroups;
            public int DistinctPeptides;
            public int DistinctMatches;

            #region Constructor
            public TotalCounts (NHibernate.ISession session, DataFilter dataFilter)
            {
                lock (session)
                {
                    var total = session.CreateQuery("SELECT " +
                                                    "COUNT(DISTINCT psm.Peptide.PeptideGroup), " +
                                                    "COUNT(DISTINCT psm.Peptide.id), " +
                                                    "COUNT(DISTINCT psm.DistinctMatchKey) " +
                                                    dataFilter.GetFilteredQueryString(DataFilter.FromPeptideSpectrumMatch))
                        .UniqueResult<object[]>();

                    PeptideGroups = Convert.ToInt32(total[0]);
                    DistinctPeptides = Convert.ToInt32(total[1]);
                    DistinctMatches = Convert.ToInt32(total[2]);
                }
            }
            #endregion
        }
        #endregion

        #region Functions for getting rows
        IList<Row> getPeptideGroupRows (DataFilter parentFilter)
        {
            lock (session)
            return session.CreateQuery(AggregateRow.Selection + ", pep.PeptideGroup " +
                                       parentFilter.GetFilteredQueryString(DataFilter.FromPeptide,
                                                                           DataFilter.PeptideToPeptideSpectrumMatch,
                                                                           DataFilter.PeptideToProtein) +
                                       "GROUP BY pep.PeptideGroup " +
                                       "ORDER BY COUNT(DISTINCT pep.id) DESC")//, COUNT(DISTINCT psm.id) DESC, COUNT(DISTINCT psm.Spectrum.id) DESC")
                          .List<object[]>()
                          .Select(o => new PeptideGroupRow(o, parentFilter) as Row)
                          .ToList();
        }

        IList<Row> getDistinctPeptideRows (DataFilter parentFilter)
        {
            lock (session)
            return session.CreateQuery(AggregateRow.Selection + ", psm.Peptide " +
                                       parentFilter.GetFilteredQueryString(DataFilter.FromPeptideSpectrumMatch, DataFilter.PeptideSpectrumMatchToProtein) +
                                       "GROUP BY psm.Peptide.Id " +
                                       "ORDER BY COUNT(DISTINCT psm.Peptide.id) DESC")//, COUNT(DISTINCT psm.id) DESC, COUNT(DISTINCT psm.Spectrum.id) DESC")
                          .List<object[]>()
                          .Select(o => new DistinctPeptideRow(o, parentFilter) as Row)
                          .ToList();
        }

        IList<Row> getDistinctMatchRows (DataFilter parentFilter)
        {
            lock (session)
            return session.CreateQuery(AggregateRow.Selection + ", psm.Peptide, psm, psm.DistinctMatchKey " +
                                       parentFilter.GetFilteredQueryString(DataFilter.FromPeptideSpectrumMatch, DataFilter.PeptideSpectrumMatchToProtein) +
                                       "GROUP BY psm.DistinctMatchKey " +
                                       "ORDER BY COUNT(DISTINCT psm.id) DESC")//, COUNT(DISTINCT psm.id) DESC, COUNT(DISTINCT psm.Spectrum.id) DESC")
                          .List<object[]>()
                          .Select(o => new DistinctMatchRow(o, parentFilter) as Row)
                          .ToList();
        }

        IList<Row> getChildren (Grouping<GroupBy> grouping, DataFilter parentFilter)
        {
            if (grouping == null)
                return getDistinctMatchRows(parentFilter);

            switch (grouping.Mode)
            {
                case GroupBy.PeptideGroup: return getPeptideGroupRows(parentFilter);
                case GroupBy.Peptide: return getDistinctPeptideRows(parentFilter);
                default: throw new NotImplementedException();
            }
        }

        protected override IList<Row> getChildren (Row parentRow)
        {
            if (parentRow.ChildRows != null)
                return parentRow.ChildRows;

            if (parentRow is PeptideGroupRow)
            {
                var row = parentRow as PeptideGroupRow;
                var parentFilter = row.DataFilter ?? dataFilter;
                var childFilter = new DataFilter(parentFilter) { PeptideGroup = new List<int>() { row.PeptideGroup } };
                var childGrouping = GroupingSetupControl<GroupBy>.GetChildGrouping(checkedGroupings, GroupBy.PeptideGroup);
                parentRow.ChildRows = getChildren(childGrouping, childFilter);
            }
            else if (parentRow is DistinctPeptideRow)
            {
                var row = parentRow as DistinctPeptideRow;
                var parentFilter = row.DataFilter ?? dataFilter;
                var childFilter = new DataFilter(parentFilter) { Peptide = new List<Peptide>() { row.Peptide } };
                var childGrouping = GroupingSetupControl<GroupBy>.GetChildGrouping(checkedGroupings, GroupBy.Peptide);
                parentRow.ChildRows = getChildren(childGrouping, childFilter);
            }
            else if (parentRow is AggregateRow)
                throw new NotImplementedException();
            else if (parentRow == null)
            {
                return getDistinctMatchRows(dataFilter);
            }

            return parentRow.ChildRows;
        }

        static string pivotHqlFormat = @"SELECT {0}, {1}, {2},
                                                COUNT(DISTINCT psm.Spectrum.id),
                                                COUNT(DISTINCT psm.DistinctMatchKey),
                                                COUNT(DISTINCT psm.Peptide.id)
                                         {3}
                                         GROUP BY {0}, {1}
                                         ORDER BY {0}
                                        ";
        Map<long, Map<long, PivotData>> getPivotData (Grouping<GroupBy> group, Pivot<PivotBy> pivot, DataFilter parentFilter)
        {
            // ProteinGroup and Cluster are consecutive, 1-based series
            string groupColumn = "psm.DistinctMatchKey";
            if (group != null)
            {
                if (group.Mode == GroupBy.PeptideGroup) groupColumn = "pep.PeptideGroup";
                else if (group.Mode == GroupBy.Peptide) groupColumn = "psm.Peptide.id";
                else throw new ArgumentException();
            }
            var pivotColumn = pivot.Text.Contains("Group") ? "ssgl.Group.id" : "s.Source.id";
            var rowIdColumn = groupColumn == "psm.DistinctMatchKey" ? "psm.id" : groupColumn;

            var pivotHql = String.Format(pivotHqlFormat,
                                         groupColumn, pivotColumn, rowIdColumn,
                                         parentFilter.GetFilteredQueryString(DataFilter.FromPeptideSpectrumMatch,
                                                                             DataFilter.PeptideSpectrumMatchToSpectrumSourceGroupLink,
                                                                             DataFilter.PeptideSpectrumMatchToPeptide));
            var query = session.CreateQuery(pivotHql);
            var pivotData = new Map<long, Map<long, PivotData>>();

            IList<object[]> pivotRows; lock (session) pivotRows = query.List<object[]>();
            foreach (var queryRow in pivotRows)
                pivotData[Convert.ToInt64(queryRow[1])][Convert.ToInt64(queryRow[2])] = new PivotData(queryRow);
            return pivotData;
        }

        #endregion

        public event PeptideViewFilterEventHandler PeptideViewFilter;

        private TotalCounts totalCounts, basicTotalCounts;

        private Dictionary<long, SpectrumSource> sourceById;
        private Dictionary<long, SpectrumSourceGroup> groupById;

        // map source/group id to row index to pivot data
        private Map<long, Map<long, PivotData>> statsBySpectrumSource, basicStatsBySpectrumSource;
        private Map<long, Map<long, PivotData>> statsBySpectrumSourceGroup, basicStatsBySpectrumSourceGroup;

        // TODO: support multiple selected objects
        string[] oldSelectionPath = new string[] { };

        public PeptideTableForm ()
        {
            InitializeComponent();

            Text = TabText = "Peptide View";

            SetDefaults();

            groupingSetupControl.GroupingChanging += groupingSetupControl_GroupingChanging;
            pivotSetupControl.PivotChanged += pivotSetupControl_PivotChanged;

            treeDataGridView.CellValueNeeded += treeDataGridView_CellValueNeeded;
            treeDataGridView.CellFormatting += treeDataGridView_CellFormatting;
            treeDataGridView.CellMouseClick += treeDataGridView_CellMouseClick;
            treeDataGridView.CellContentClick += treeDataGridView_CellContentClick;
            treeDataGridView.CellDoubleClick += treeDataGridView_CellDoubleClick;
            treeDataGridView.PreviewKeyDown += treeDataGridView_PreviewKeyDown;
            treeDataGridView.CellIconNeeded += treeDataGridView_CellIconNeeded;
        }

        void treeDataGridView_CellIconNeeded (object sender, TreeDataGridViewCellValueEventArgs e)
        {
            if (e.RowIndexHierarchy.First() >= rows.Count)
            {
                e.Value = null;
                return;
            }

            Row baseRow = rows[e.RowIndexHierarchy.First()];
            for (int i = 1; i < e.RowIndexHierarchy.Count; ++i)
            {
                getChildren(baseRow); // populate ChildRows if necessary
                baseRow = baseRow.ChildRows[e.RowIndexHierarchy[i]];
            }

            if (baseRow is PeptideGroupRow) e.Value = Properties.Resources.PeptideGroup;
            else if (baseRow is DistinctPeptideRow) e.Value = Properties.Resources.Peptide;
            else if (baseRow is DistinctMatchRow) e.Value = Properties.Resources.DistinctMatch;
        }

        private void treeDataGridView_CellValueNeeded (object sender, TreeDataGridViewCellValueEventArgs e)
        {
            var rootGrouping = checkedGroupings.Count > 0 ? checkedGroupings.First() : null;

            if (e.RowIndexHierarchy.First() >= rows.Count)
            {
                e.Value = null;
                return;
            }

            Row baseRow = rows[e.RowIndexHierarchy.First()];
            for (int i = 1; i < e.RowIndexHierarchy.Count; ++i)
            {
                getChildren(baseRow); // populate ChildRows if necessary
                baseRow = baseRow.ChildRows[e.RowIndexHierarchy[i]];
            }

            if (baseRow is PeptideGroupRow)
            {
                var row = baseRow as PeptideGroupRow;
                
                var childGrouping = GroupingSetupControl<GroupBy>.GetChildGrouping(checkedGroupings, GroupBy.PeptideGroup);
                if (childGrouping == null)
                    e.ChildRowCount = row.DistinctMatches;
                else if (childGrouping.Mode == GroupBy.Peptide)
                    e.ChildRowCount = row.DistinctPeptides;
            }
            else if (baseRow is DistinctPeptideRow)
            {
                var row = baseRow as DistinctPeptideRow;
                e.ChildRowCount = row.DistinctMatches;
            }

            e.Value = getCellValue(e.ColumnIndex, baseRow);
        }

        protected override object getCellValue (int columnIndex, Row baseRow)
        {
            // TODO: fix child rows so they have their own pivot data
            if (pivotColumns.Count > 0 && columnIndex >= pivotColumns.First().Index)
            {
                var stats = treeDataGridView.Columns[columnIndex].Tag as Map<long, PivotData>;
                if (stats == null)
                    return 0;

                long rowId;
                if (baseRow is PeptideGroupRow)
                    rowId = (baseRow as PeptideGroupRow).PeptideGroup;
                else if (baseRow is DistinctPeptideRow)
                    rowId = (baseRow as DistinctPeptideRow).Peptide.Id.Value;
                else if (baseRow is DistinctMatchRow)
                    rowId = (baseRow as DistinctMatchRow).PeptideSpectrumMatch.Id.Value;
                else
                    throw new NotImplementedException();
                var itr = stats.Find(rowId);

                if (itr.IsValid)
                {
                    if (treeDataGridView.Columns[columnIndex].HeaderText.EndsWith("/"))
                    {
                        if (checkedPivots.Count(o => o.Mode == PivotBy.SpectraByGroup) > 0)
                            return itr.Current.Value.Spectra;
                        else if (checkedPivots.Count(o => o.Mode == PivotBy.MatchesByGroup) > 0)
                            return itr.Current.Value.DistinctMatches;
                        else if (checkedPivots.Count(o => o.Mode == PivotBy.PeptidesByGroup) > 0)
                            return itr.Current.Value.DistinctPeptides;
                    }
                    else
                    {
                        if (checkedPivots.Count(o => o.Mode == PivotBy.SpectraBySource) > 0)
                            return itr.Current.Value.Spectra;
                        else if (checkedPivots.Count(o => o.Mode == PivotBy.MatchesBySource) > 0)
                            return itr.Current.Value.DistinctMatches;
                        else if (checkedPivots.Count(o => o.Mode == PivotBy.PeptidesBySource) > 0)
                            return itr.Current.Value.DistinctPeptides;
                    }
                }
            }
            else if (baseRow is PeptideGroupRow)
            {
                var row = baseRow as PeptideGroupRow;
                if (columnIndex == keyColumn.Index) return row.PeptideGroup;
                else if (columnIndex == distinctPeptidesColumn.Index) return row.DistinctPeptides;
                else if (columnIndex == distinctMatchesColumn.Index) return row.DistinctMatches;
                else if (columnIndex == filteredSpectraColumn.Index) return row.Spectra;
                else if (columnIndex == proteinsColumn.Index) return row.Proteins;
                else if (columnIndex == proteinAccessionsColumn.Index) return row.ProteinAccessions;
                else if (columnIndex == proteinGroupsColumn.Index) return row.ProteinGroups;
            }
            else if (baseRow is DistinctPeptideRow)
            {
                var row = baseRow as DistinctPeptideRow;
                if (columnIndex == keyColumn.Index) return row.Peptide.Sequence;
                else if (columnIndex == distinctMatchesColumn.Index) return row.DistinctMatches;
                else if (columnIndex == filteredSpectraColumn.Index) return row.Spectra;
                else if (columnIndex == monoisotopicMassColumn.Index) return row.Peptide.MonoisotopicMass;
                else if (columnIndex == molecularWeightColumn.Index) return row.Peptide.MolecularWeight;
                else if (columnIndex == peptideGroupColumn.Index) return row.Peptide.PeptideGroup;
                else if (columnIndex == proteinsColumn.Index) return row.Proteins;
                else if (columnIndex == proteinAccessionsColumn.Index) return row.ProteinAccessions;
                else if (columnIndex == proteinGroupsColumn.Index) return row.ProteinGroups;
            }
            else if (baseRow is DistinctMatchRow)
            {
                var row = baseRow as DistinctMatchRow;
                if (columnIndex == keyColumn.Index) return row.DistinctMatch;
                else if (columnIndex == monoisotopicMassColumn.Index) return row.PeptideSpectrumMatch.MonoisotopicMass;
                else if (columnIndex == molecularWeightColumn.Index) return row.PeptideSpectrumMatch.MolecularWeight;
                else if (columnIndex == filteredSpectraColumn.Index) return row.Spectra;
                else if (columnIndex == peptideGroupColumn.Index) return row.Peptide.PeptideGroup;
                else if (columnIndex == proteinsColumn.Index) return row.Proteins;
                else if (columnIndex == proteinAccessionsColumn.Index) return row.ProteinAccessions;
                else if (columnIndex == proteinGroupsColumn.Index) return row.ProteinGroups;
            }
            return null;
        }

        protected override RowFilterState getRowFilterState (Row parentRow)
        {
            bool result = false;
            if (parentRow is PeptideGroupRow)
            {
                if (viewFilter.PeptideGroup != null) result = viewFilter.PeptideGroup.Contains((parentRow as PeptideGroupRow).PeptideGroup);
            }
            else if (parentRow is DistinctPeptideRow)
            {
                if (viewFilter.Peptide != null) result = viewFilter.Peptide.Contains((parentRow as DistinctPeptideRow).Peptide);
                if (!result && viewFilter.PeptideGroup != null) result = viewFilter.PeptideGroup.Contains((parentRow as DistinctPeptideRow).Peptide.PeptideGroup);
            }
            else if (parentRow is DistinctMatchRow)
            {
                if (viewFilter.DistinctMatchKey != null) result = viewFilter.DistinctMatchKey.Contains((parentRow as DistinctMatchRow).DistinctMatch);
                if (!result && viewFilter.Peptide != null) result = viewFilter.Peptide.Contains((parentRow as DistinctMatchRow).Peptide);
                if (!result && viewFilter.PeptideGroup != null) result = viewFilter.PeptideGroup.Contains((parentRow as DistinctMatchRow).Peptide.PeptideGroup);
                result = result || viewFilter.PeptideGroup == null && viewFilter.Peptide == null && viewFilter.DistinctMatchKey == null;
            }
            if (result) return RowFilterState.In;
            if (parentRow.ChildRows == null) return RowFilterState.Out;

            return parentRow.ChildRows.Aggregate(RowFilterState.Unknown, (x, y) => x | getRowFilterState(y));
        }

        private void treeDataGridView_CellFormatting (object sender, TreeDataGridViewCellFormattingEventArgs e)
        {
            var column = treeDataGridView.Columns[e.ColumnIndex];
            if (_columnSettings.ContainsKey(column) && _columnSettings[column].BackColor.HasValue)
                e.CellStyle.BackColor = _columnSettings[column].BackColor.Value;
            else
                e.CellStyle.BackColor = e.CellStyle.BackColor;

            if (viewFilter.Cluster == null && viewFilter.Peptide == null && viewFilter.DistinctMatchKey == null)
                return;

            Row row = rows[e.RowIndexHierarchy.First()];
            for (int i = 1; i < e.RowIndexHierarchy.Count; ++i)
                row = row.ChildRows[e.RowIndexHierarchy[i]];

            switch (getRowFilterState(row))
            {
                case RowFilterState.Out:
                    e.CellStyle.ForeColor = filteredOutColor;
                    break;
                case RowFilterState.Partial:
                    e.CellStyle.ForeColor = filteredPartialColor;
                    break;
            }

            if (column is DataGridViewLinkColumn)
            {
                var cell = treeDataGridView[e.ColumnIndex, e.RowIndexHierarchy] as DataGridViewLinkCell;
                cell.LinkColor = cell.ActiveLinkColor = e.CellStyle.ForeColor;
            }
        }

        void treeDataGridView_CellMouseClick (object sender, TreeDataGridViewCellMouseEventArgs e)
        {
            if (e.ColumnIndex < 0)
                return;

            // was column header clicked?
            if (e.RowIndexHierarchy.First() < 0)
                Sort(e.ColumnIndex);
        }

        void treeDataGridView_CellContentClick (object sender, TreeDataGridViewCellEventArgs e)
        {
            if (e.ColumnIndex < 0 || e.RowIndexHierarchy.First() < 0)
                return;

            /*if (e.ColumnIndex == clusterColumn.Index && ProteinViewFilter != null)
            {
                object value = treeDataGridView[e.ColumnIndex, e.RowIndexHierarchy].Value;
                if (value == null)
                    return;

                var newDataFilter = new DataFilter(dataFilter) { FilterSource = this };
                newDataFilter.Cluster = new List<int> { (int) value };

                PeptideViewFilter(this, newDataFilter);
            }
            else if (e.ColumnIndex == coverageColumn.Index && ProteinViewVisualize != null)
            {
                Row row = rows[e.RowIndexHierarchy.First()];
                for (int i = 1; i < e.RowIndexHierarchy.Count; ++i)
                    row = row.ChildRows[e.RowIndexHierarchy[i]];

                if (row is ProteinRow)
                    ProteinViewVisualize(this, new ProteinViewVisualizeEventArgs() { Protein = (row as ProteinRow).Protein });
            }*/
        }

        void treeDataGridView_CellDoubleClick (object sender, TreeDataGridViewCellEventArgs e)
        {
            if (e.ColumnIndex < 0 || e.RowIndexHierarchy.First() < 0)
                return;

            Row row = rows[e.RowIndexHierarchy.First()];
            for (int i = 1; i < e.RowIndexHierarchy.Count; ++i)
                row = row.ChildRows[e.RowIndexHierarchy[i]];

            var newDataFilter = new DataFilter() { FilterSource = this };

            if (row is PeptideGroupRow)
                newDataFilter.PeptideGroup = new List<int>() { (row as PeptideGroupRow).PeptideGroup };
            if (row is DistinctPeptideRow)
                newDataFilter.Peptide = new List<Peptide>() { (row as DistinctPeptideRow).Peptide };
            else if (row is DistinctMatchRow)
                newDataFilter.DistinctMatchKey = new List<DistinctMatchKey>() { (row as DistinctMatchRow).DistinctMatch };

            if (PeptideViewFilter != null)
                PeptideViewFilter(this, newDataFilter);
        }

        void treeDataGridView_PreviewKeyDown (object sender, PreviewKeyDownEventArgs e)
        {
            if (e.KeyCode == Keys.Escape)
                treeDataGridView.ClearSelection();

            if (e.KeyCode != Keys.Enter)
                return;

            var newDataFilter = new DataFilter { FilterSource = this };

            if (treeDataGridView.SelectedCells.Count == 0)
                return;

            var processedRows = new Set<int>();
            var selectedPeptideGroups = new List<int>();
            var selectedPeptides = new List<Peptide>();
            var selectedMatches = new List<DistinctMatchKey>();

            foreach (DataGridViewCell cell in treeDataGridView.SelectedCells)
            {
                if (!processedRows.Insert(cell.RowIndex).WasInserted)
                    continue;

                var rowIndexHierarchy = treeDataGridView.GetRowHierarchyForRowIndex(cell.RowIndex);

                Row row = rows[rowIndexHierarchy.First()];
                for (int i = 1; i < rowIndexHierarchy.Count; ++i)
                    row = row.ChildRows[rowIndexHierarchy[i]];

                if (row is PeptideGroupRow)
                    selectedPeptideGroups.Add((row as PeptideGroupRow).PeptideGroup);
                else if (row is DistinctPeptideRow)
                    selectedPeptides.Add((row as DistinctPeptideRow).Peptide);
                else if (row is DistinctMatchRow)
                    selectedMatches.Add((row as DistinctMatchRow).DistinctMatch);
            }

            if (selectedPeptideGroups.Count > 0) newDataFilter.PeptideGroup = selectedPeptideGroups;
            if (selectedPeptides.Count > 0) newDataFilter.Peptide = selectedPeptides;
            if (selectedMatches.Count > 0) newDataFilter.DistinctMatchKey = selectedMatches;

            if (PeptideViewFilter != null)
                PeptideViewFilter(this, newDataFilter);
        }

        private void SetDefaults()
        {
            _columnSettings = new Dictionary<DataGridViewColumn, ColumnProperty>()
            {
                { keyColumn, new ColumnProperty() {Type = typeof(string)}},
                { distinctPeptidesColumn, new ColumnProperty() {Type = typeof(int)}},
                { distinctMatchesColumn, new ColumnProperty() {Type = typeof(int)}},
                { filteredSpectraColumn, new ColumnProperty() {Type = typeof(int)}},
                { monoisotopicMassColumn, new ColumnProperty() {Type = typeof(float)}},
                { molecularWeightColumn, new ColumnProperty() {Type = typeof(float)}},
                { peptideGroupColumn, new ColumnProperty() {Type = typeof(int), Visible = false}},
                { proteinsColumn, new ColumnProperty() {Type = typeof(int), Visible = false}},
                { proteinAccessionsColumn, new ColumnProperty() {Type = typeof(string), Visible = false}},
                { proteinGroupsColumn, new ColumnProperty() {Type = typeof(string), Visible = false}},
            };

            foreach (var kvp in _columnSettings)
            {
                kvp.Value.Name = kvp.Key.Name;
                kvp.Value.Index = kvp.Key.Index;
                kvp.Value.DisplayIndex = kvp.Key.DisplayIndex;
            }

            initialColumnSortOrders = new Map<int, SortOrder>()
            {
                {keyColumn.Index, SortOrder.Ascending},
                {distinctPeptidesColumn.Index, SortOrder.Descending},
                {distinctMatchesColumn.Index, SortOrder.Descending},
                {filteredSpectraColumn.Index, SortOrder.Descending},
                {monoisotopicMassColumn.Index, SortOrder.Ascending},
                {molecularWeightColumn.Index, SortOrder.Ascending},
                {peptideGroupColumn.Index, SortOrder.Ascending},
                {proteinsColumn.Index, SortOrder.Ascending},
                {proteinAccessionsColumn.Index, SortOrder.Ascending},
                {proteinGroupsColumn.Index, SortOrder.Ascending},
            };
        }

        public override void SetData (NHibernate.ISession session, DataFilter dataFilter)
        {
            this.session = session;
            viewFilter = dataFilter;
            this.dataFilter = new DataFilter(dataFilter) {Peptide = null, DistinctMatchKey = null};

            /*if (treeListView.SelectedObject is PeptideRow)
                oldSelectionPath = new string[] { treeListView.SelectedItem.Text };
            else if (treeListView.SelectedObject is PeptideSpectrumMatchRow)
                oldSelectionPath = new string[] { (treeListView.SelectedObject as PeptideSpectrumMatchRow).PeptideSpectrumMatch.Peptide.Sequence, treeListView.SelectedItem.Text };*/

            ClearData();

            // stored to avoid cross-thread calls on the control
            checkedPivots = pivotSetupControl.CheckedPivots;
            checkedGroupings = groupingSetupControl.CheckedGroupings;

            Text = TabText = "Loading peptide view...";

            var workerThread = new BackgroundWorker()
            {
                WorkerReportsProgress = true,
                WorkerSupportsCancellation = true
            };

            workerThread.DoWork += new DoWorkEventHandler(setData);
            workerThread.RunWorkerCompleted += new RunWorkerCompletedEventHandler(renderData);
            workerThread.RunWorkerAsync();
        }

        private void addPivotColumns ()
        {
            treeDataGridView.SuspendLayout();
            foreach (var pivotColumn in pivotColumns)
                treeDataGridView.Columns.Remove(pivotColumn);

            pivotColumns = new List<DataGridViewColumn>();

            if (checkedPivots.Count == 0)
            {
                treeDataGridView.ResumeLayout(true);
                return;
            }

            var sourceNames = sourceById.Select(o => o.Value.Name);
            int insertIndex = treeDataGridView.GetVisibleColumns().Count();

            if (statsBySpectrumSource != null)
                foreach (long sourceId in statsBySpectrumSource.Keys)
                {
                    string uniqueSubstring;
                    Util.UniqueSubstring(sourceById[sourceId].Name, sourceNames, out uniqueSubstring);
                    var column = new DataGridViewTextBoxColumn() { HeaderText = uniqueSubstring };
                    column.Tag = statsBySpectrumSource[sourceId];

                    var newProperties = new ColumnProperty()
                    {
                        Type = typeof(int),
                        Name = column.HeaderText
                    };

                    var previousForm = _columnSettings.SingleOrDefault(x => x.Value.Name == column.HeaderText);

                    if (previousForm.Key != null)
                    {
                        _columnSettings.Remove(previousForm.Key);
                        newProperties = previousForm.Value;
                    }
                    else
                    {
                        var possibleSaved = _unusedPivotSettings.SingleOrDefault(x => x.Name == column.HeaderText);
                        if (possibleSaved != null)
                            newProperties = possibleSaved;
                    }

                    if (newProperties.Visible.HasValue)
                        column.Visible = newProperties.Visible.Value;
                    //_columnSettings.Add(column, newProperties);
                    if (column.Visible)
                        pivotColumns.Add(column);
                }

            if (statsBySpectrumSourceGroup != null)
                foreach (long groupId in statsBySpectrumSourceGroup.Keys)
                {
                    var column = new DataGridViewTextBoxColumn() { HeaderText = groupById[groupId].Name.TrimEnd('/') + '/' };
                    column.Tag = statsBySpectrumSourceGroup[groupId];

                    var newProperties = new ColumnProperty()
                    {
                        Type = typeof(int),
                        Name = column.HeaderText
                    };

                    var previousForm = _columnSettings.SingleOrDefault(x => x.Value.Name == column.HeaderText);

                    if (previousForm.Key != null)
                    {
                        _columnSettings.Remove(previousForm.Key);
                        newProperties = previousForm.Value;
                    }
                    else
                    {
                        var possibleSaved = _unusedPivotSettings.SingleOrDefault(x => x.Name == column.HeaderText);
                        if (possibleSaved != null)
                            newProperties = possibleSaved;
                    }

                    if (newProperties.Visible.HasValue)
                        column.Visible = newProperties.Visible.Value;
                    //_columnSettings.Add(column, newProperties);
                    if (column.Visible)
                        pivotColumns.Add(column);
                }

            pivotColumns.Sort((x, y) => x.HeaderText.CompareTo(y.HeaderText));
            treeDataGridView.Columns.AddRange(pivotColumns.ToArray());
            treeDataGridView.ResumeLayout(true);
        }

        List<ColumnProperty> _unusedPivotSettings = new List<ColumnProperty>();

        protected override bool updatePivots (FormProperty formProperty)
        {
            if (pivotSetupControl != null && formProperty.PivotModes != null)
                return base.updatePivots(formProperty);
            else
            {
                setPivots(new Pivot<PivotBy>() {Mode = PivotBy.SpectraByGroup, Text = "Spectra by Group"},
                          new Pivot<PivotBy>() {Mode = PivotBy.SpectraBySource, Text = "Spectra by Source"},
                          new Pivot<PivotBy>() {Mode = PivotBy.MatchesByGroup, Text = "Matches by Group"},
                          new Pivot<PivotBy>() {Mode = PivotBy.MatchesBySource, Text = "Matches by Source"},
                          new Pivot<PivotBy>() {Mode = PivotBy.PeptidesByGroup, Text = "Peptides by Group"},
                          new Pivot<PivotBy>() {Mode = PivotBy.PeptidesBySource, Text = "Peptides by Source"});
                pivotSetupControl.PivotChanged += pivotSetupControl_PivotChanged;
                return false;
            }
        }

        protected override bool updateGroupings (FormProperty formProperty)
        {
            bool groupingChanged = false;
            if (groupingSetupControl != null && formProperty.GroupingModes != null)
                groupingChanged = base.updateGroupings(formProperty);
            else
                setGroupings(new Grouping<GroupBy>() { Mode = GroupBy.PeptideGroup, Text = "Peptide Group" },
                             new Grouping<GroupBy>() { Mode = GroupBy.Peptide, Text = "Peptide" });

            groupingSetupControl.GroupingChanging += groupingSetupControl_GroupingChanging;

            if (groupingChanged)
                setColumnVisibility();

            return groupingChanged;
        }

        public override void ClearData ()
        {
            Text = TabText = "Peptide View";

            treeDataGridView.RootRowCount = 0;
            Refresh();
        }

        public override void ClearData (bool clearBasicFilter)
        {
            if (clearBasicFilter)
                basicDataFilter = null;
            ClearData();
        }

        void setData(object sender, DoWorkEventArgs e)
        {
            try
            {
                var rootGrouping = checkedGroupings.Count > 0 ? checkedGroupings.First() : null;

                if (dataFilter.IsBasicFilter)
                {
                    // refresh basic data when basicDataFilter is unset or when the basic filter values have changed
                    if (basicDataFilter == null || (viewFilter.IsBasicFilter && dataFilter != basicDataFilter))
                    {
                        basicDataFilter = dataFilter;
                        basicTotalCounts = new TotalCounts(session, dataFilter);
                        basicRows = getChildren(rootGrouping, dataFilter);

                        basicStatsBySpectrumSourceGroup = null;
                        Pivot<PivotBy> pivotBySource = checkedPivots.FirstOrDefault(o => o.Mode.ToString().Contains("Source"));
                        if (pivotBySource != null)
                            basicStatsBySpectrumSource = getPivotData(rootGrouping, pivotBySource, dataFilter);

                        basicStatsBySpectrumSourceGroup = null;
                        Pivot<PivotBy> pivotByGroup = checkedPivots.FirstOrDefault(o => o.Mode.ToString().Contains("Group"));
                        if (pivotByGroup != null)
                            basicStatsBySpectrumSourceGroup = getPivotData(rootGrouping, pivotByGroup, dataFilter);

                        lock (session)
                        {
                            sourceById = session.Query<SpectrumSource>().Where(o => o.Group != null).ToDictionary(o => o.Id.Value);
                            groupById = session.Query<SpectrumSourceGroup>().ToDictionary(o => o.Id.Value);
                        }
                    }

                    totalCounts = basicTotalCounts;
                    rows = basicRows;
                    statsBySpectrumSource = basicStatsBySpectrumSource;
                    statsBySpectrumSourceGroup = basicStatsBySpectrumSourceGroup;
                }
                else
                {
                    totalCounts = new TotalCounts(session, dataFilter);
                    rows = getChildren(rootGrouping, dataFilter);

                    statsBySpectrumSource = null;
                    Pivot<PivotBy> pivotBySource = checkedPivots.FirstOrDefault(o => o.Mode.ToString().Contains("Source"));
                    if (pivotBySource != null)
                        statsBySpectrumSource = getPivotData(rootGrouping, pivotBySource, dataFilter);

                    statsBySpectrumSourceGroup = null;
                    Pivot<PivotBy> pivotByGroup = checkedPivots.FirstOrDefault(o => o.Mode.ToString().Contains("Group"));
                    if (pivotByGroup != null)
                        statsBySpectrumSourceGroup = getPivotData(rootGrouping, pivotByGroup, dataFilter);
                }

                applySort();
            }
            catch (Exception ex)
            {
                e.Result = ex;
            }
        }

        void renderData (object sender, RunWorkerCompletedEventArgs e)
        {
            if (e.Result is Exception)
                Program.HandleException(e.Result as Exception);

            treeDataGridView.RootRowCount = rows.Count();

            // show total counts in the form title
            Text = TabText = String.Format("Peptide View: {0} peptide groups, {1} distinct peptides, {2} distinct matches",
                                           totalCounts.PeptideGroups, totalCounts.DistinctPeptides, totalCounts.DistinctMatches);

            addPivotColumns();

            // try to (re)set selected item
            expandSelectionPath(oldSelectionPath);

            treeDataGridView.Refresh();
        }

        private void expandSelectionPath (IEnumerable<string> selectionPath)
        {
            /*OLVListItem selectedItem = null;
            foreach (string branch in selectionPath)
            {
                int index = 0;
                if (selectedItem != null)
                {
                    treeListView.Expand(selectedItem.RowObject);
                    index = selectedItem.Index;
                }

                index = treeListView.FindMatchingRow(branch, index, SearchDirectionHint.Down);
                if (index < 0)
                    break;
                selectedItem = treeListView.Items[index] as OLVListItem;
            }

            if (selectedItem != null)
            {
                treeListView.SelectedItem = selectedItem;
                selectedItem.EnsureVisible();
            }*/
        }

        private void groupingSetupControl_GroupingChanging (object sender, GroupingChangingEventArgs<GroupBy> e)
        {
            // GroupBy.Peptide cannot be before GroupBy.PeptideGroup

            if (e.Grouping.Mode != GroupBy.PeptideGroup && e.Grouping.Mode != GroupBy.Peptide)
                return;

            var newGroupings = new List<Grouping<GroupBy>>(groupingSetupControl.Groupings);
            newGroupings.Remove(e.Grouping);
            newGroupings.Insert(e.NewIndex, e.Grouping);

            e.Cancel = GroupingSetupControl<GroupBy>.HasParentGrouping(newGroupings, GroupBy.PeptideGroup, GroupBy.Peptide);
        }

        protected override void setColumnVisibility ()
        {
            var columnsIrrelevantForGrouping = new Set<DataGridViewColumn>(new Comparison<DataGridViewColumn>((x, y) => x.Name.CompareTo(y.Name)));
            if (checkedGroupings.IsNullOrEmpty())
            {
                columnsIrrelevantForGrouping.Add(distinctPeptidesColumn);
                columnsIrrelevantForGrouping.Add(distinctMatchesColumn);
            }
            else if (checkedGroupings.First().Mode == GroupBy.PeptideGroup)
                columnsIrrelevantForGrouping.Add(peptideGroupColumn);
            else if (checkedGroupings.First().Mode == GroupBy.Peptide)
                columnsIrrelevantForGrouping.Add(distinctPeptidesColumn);

            // if visibility is not forced, use grouping mode to set automatic visibility
            foreach (var kvp in _columnSettings)
                kvp.Key.Visible = kvp.Value.Visible ?? !columnsIrrelevantForGrouping.Contains(kvp.Key);

            base.setColumnVisibility();
        }

        protected override void OnGroupingChanged (object sender, EventArgs e)
        {
            setColumnVisibility();
            base.OnGroupingChanged(sender, e);
        }

        private void pivotSetupControl_PivotChanged (object sender, PivotChangedEventArgs<PivotBy> e)
        {
            if (e.Pivot.Checked)
            {
                // uncheck mutually exclusive pivot modes
                string exclusiveMode = e.Pivot.Text.Contains("Source") ? "Source" : "Group";
                var conflictingPivots = pivotSetupControl.CheckedPivots.Where(o => o != e.Pivot && o.Text.Contains(exclusiveMode));
                foreach (var pivot in conflictingPivots)
                    pivotSetupControl.SetPivot(pivot.Mode, false);
            }
        }
    }

    public delegate void PeptideViewFilterEventHandler (PeptideTableForm sender, DataFilter peptideViewFilter);
}