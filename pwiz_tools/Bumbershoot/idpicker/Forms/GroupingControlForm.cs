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
using System.Collections;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Windows.Forms;
using BrightIdeasSoftware;
using NHibernate.Linq;
using IDPicker.DataModel;

namespace IDPicker.Forms
{
    public partial class GroupingControlForm : Form
    {
        private class tlvBranch
        {
            public string Text;
            public object Data;
            public tlvBranch Parent;
            public List<object> Children;
            public ContextMenuStrip cms;
        }

        NHibernate.ISession session;
        private OLVListItem _clickedItem;
        int _numberNewNodes;
        private tlvBranch _rootNode;

        public GroupingControlForm(NHibernate.ISessionFactory sessionFactory)
        {
            InitializeComponent();

            this.session = sessionFactory.OpenSession();

            tlvGroupedFiles.CanExpandGetter += x => (((tlvBranch)x).Data is SpectrumSourceGroup && ((tlvBranch)x).Children.Any());
            tlvGroupedFiles.ChildrenGetter += getChildren;
            tlvGroups.AspectGetter += x =>
                                          {
                                              var nodeTracker = ((tlvBranch) x);
                                              var offsetCorrection = 0;
                                              while (nodeTracker.Parent.Text != null)
                                              {
                                                  offsetCorrection++;
                                                  nodeTracker = nodeTracker.Parent;
                                              }
                                              return ((tlvBranch) x).Text + new string(' ',offsetCorrection*7);
                                          };
            tlvGroups.ImageGetter += delegate(object x) { return (((tlvBranch)x).Data is SpectrumSourceGroup) ? Properties.Resources.XPfolder_closed : Properties.Resources.file; };

            ApplyDefaultGroups(null, null);
        }

        private IEnumerable getChildren(object model)
        {
            var branch = (tlvBranch)model;
            return branch.Children;
        }

        private void AddAssociatedSpectra(ref tlvBranch groupNode, Dictionary<int, List<SpectrumSource>> spectraDictionary)
        {
            try
            {
                int groupID = (int)(groupNode.Data as SpectrumSourceGroup).Id;

                if (spectraDictionary.ContainsKey(groupID))
                {
                    foreach (SpectrumSource ss in spectraDictionary[groupID])
                    {
                        var newNode = new tlvBranch
                                          {
                                              Text = Path.GetFileName(ss.Name),
                                              Parent = groupNode,
                                              Data = ss,
                                              cms = cmRightClickFileNode
                                          };

                        groupNode.Children.Add(newNode);
                    }
                }
            }
            catch
            {
                MessageBox.Show("Could not add spectra to branch");
            }
        }

        private tlvBranch FillBranch(tlvBranch groupNode, IList<SpectrumSourceGroup> groups, Dictionary<int, List<SpectrumSource>> spectraDictionary)
        {
            List<SpectrumSourceGroup> potentialChildren;
            string fullPath = (groupNode.Data as SpectrumSourceGroup).Name;

            if (fullPath == "/")
                potentialChildren = (from g in groups
                                     where g.Name != fullPath &&
                                     !g.Name.Remove(0,1).Contains("/")
                                     select g).ToList();
            else
                potentialChildren = (from g in groups
                                     where g.Name.Contains(fullPath + "/") &&
                                     !g.Name.Remove(0,fullPath.Length+1).Contains("/")
                                     select g).ToList();

            foreach (SpectrumSourceGroup ssg in potentialChildren)
            {
                var newNode = new tlvBranch
                                  {
                                      Text = Path.GetFileName(ssg.Name),
                                      Children = new List<object>(),
                                      Data = ssg,
                                      Parent = groupNode,
                                      cms = cmRightClickGroupNode
                                  };
                groups.Remove(newNode.Data as SpectrumSourceGroup);

                AddAssociatedSpectra(ref newNode, spectraDictionary);

                groupNode.Children.Add(FillBranch(newNode, groups, spectraDictionary));
            }

            OrganizeNode(groupNode);
            return groupNode;
        }

        private void saveButton_Click(object sender, EventArgs e)
        {
            var spectraLocations = new List<tlvBranch>();
            var groupsToSave = new List<tlvBranch>();

            var transaction = session.BeginTransaction();

            // find all groups still present
            getSprectrumSourceGroupsRecursively(_rootNode, groupsToSave,string.Empty);

            // remove old groups and links
            session.CreateQuery("DELETE SpectrumSourceGroupLink").ExecuteUpdate();
            var unusedGroups = session.Query<SpectrumSourceGroup>().ToList();

            foreach (tlvBranch tn in groupsToSave)
                unusedGroups.Remove((tn.Data as SpectrumSourceGroup));
            foreach (SpectrumSourceGroup ssg in unusedGroups)
                session.Delete(ssg);
               
            // save group layout
            foreach (tlvBranch treeNode in groupsToSave)
                session.SaveOrUpdate(treeNode.Data as SpectrumSourceGroup);

            // get new spectra locations
            getListOfSprectrumSourcesRecursively(_rootNode, ref spectraLocations);

            // update SpectrumSource.Group_ and insert new SpectrumSourceGroupLinks;
            // using prepared SQL commands for speed

            var cmd1 = session.Connection.CreateCommand();
            var cmd2 = session.Connection.CreateCommand();
            cmd1.CommandText = "UPDATE SpectrumSource SET Group_ = ? WHERE Id = ?";
            cmd2.CommandText = "INSERT INTO SpectrumSourceGroupLink (Group_, Source) VALUES (?,?)";
            var parameters1 = new List<IDbDataParameter>();
            var parameters2 = new List<IDbDataParameter>();
            for (int i = 0; i < 2; ++i)
            {
                parameters1.Add(cmd1.CreateParameter());
                parameters2.Add(cmd2.CreateParameter());
                cmd1.Parameters.Add(parameters1[i]);
                cmd2.Parameters.Add(parameters2[i]);
            }
            cmd1.Prepare();
            cmd2.Prepare();

            foreach (tlvBranch tn in spectraLocations)
            {
                parameters1[0].Value = ((SpectrumSourceGroup)tn.Parent.Data).Id;
                parameters1[1].Value = ((SpectrumSource)tn.Data).Id;
                cmd1.ExecuteNonQuery();

                var tempNode = tn;
                while (tempNode.Parent.Text != null)
                {
                    parameters2[0].Value = ((SpectrumSourceGroup)tempNode.Parent.Data).Id;
                    parameters2[1].Value = ((SpectrumSource)tn.Data).Id;
                    tempNode = tempNode.Parent;
                    cmd2.ExecuteNonQuery();
                }
            }

            // save ungrouped spectrum sources
            foreach (ListViewItem lvi in lvNonGroupedFiles.Items)
            {
                var ss = (SpectrumSource)lvi.Tag;
                ss.Group = null;
                session.Update(ss);
            }

            transaction.Commit();
        }

        /// <summary>
        /// Allow delete (delete key) and rename (F2 key) of group nodes (not root or file)
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void tvGroups_KeyDown(object sender, KeyEventArgs e)
        {
            var selNode = (tlvBranch)tlvGroupedFiles.SelectedObject;

            try
            {
                if (e.KeyCode == Keys.Delete)
                    RemoveGroupNode(selNode);
            }
            catch (Exception exc)
            {
                //HandleExceptions(exc, ExceptionsDialogForm.ExceptionType.Error);
                MessageBox.Show(exc.ToString());
            }
        }

        /// <summary>
        /// Set location for determining node that was clicked on
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void tvGroups_MouseDown(object sender, MouseEventArgs e)
        {
            try
            {
                if (e.Button == MouseButtons.Right)
                {
                    var column = (OLVColumn)tlvGroupedFiles.Columns[0].Clone();
                    _clickedItem = tlvGroupedFiles.GetItemAt(e.X, e.Y, out column); ;
                    if (_clickedItem != null)
                    {
                        var selNode = (tlvBranch)tlvGroupedFiles.GetModelObject(_clickedItem.Index);
                        if (selNode != null)
                            selNode.cms.Show(tlvGroupedFiles, e.Location.X, e.Location.Y);
                    }
                }
            }
            catch (Exception exc)
            {
                //HandleExceptions(exc, ExceptionsDialogForm.ExceptionType.Error);
                MessageBox.Show(exc.ToString());
            }
        }

        /// <summary>
        /// Check to prevent group nodes from being dragged into their children
        /// group nodes.
        /// </summary>
        /// <returns></returns>
        private bool checkIfDestGroupAChildNodeOfMe(tlvBranch destNode, tlvBranch dragNode)
        {
            try
            {
                return destNode.Data is SpectrumSourceGroup
                       && dragNode.Data is SpectrumSourceGroup
                       && destNode.Text
                              .Contains(dragNode.Text);
            }
            catch (Exception exc)
            {
                throw new Exception("Error validating group move: checkIfDestGroupAChildNodeOfMe()\r\n", exc);
            }

        }

        private void getListOfSprectrumSourcesRecursively(tlvBranch treeNode, ref List<tlvBranch> spectraNodes)
        {
            try
            {
                if (treeNode.Data is SpectrumSource)
                    spectraNodes.Add(treeNode);
                else
                    foreach (tlvBranch subNode in treeNode.Children)
                        getListOfSprectrumSourcesRecursively(subNode, ref spectraNodes);
            }
            catch (Exception exc)
            {
                throw new Exception("Error in retrieving spectra", exc);
            }

        }

        private void getSprectrumSourceGroupsRecursively(tlvBranch treeNode, List<tlvBranch> nodeList, string rootname)
        {
            if (treeNode.Data is SpectrumSourceGroup)
            {
                (treeNode.Data as SpectrumSourceGroup).Name =
                    (rootname + treeNode.Text)
                        .Replace("//", "/");

                nodeList.Add(treeNode);

                foreach (tlvBranch subNode in treeNode.Children)
                    getSprectrumSourceGroupsRecursively(subNode, nodeList, (treeNode.Data as SpectrumSourceGroup).Name);
            }
        }

        private void addGroupToolStripMenuItem_Click(object sender, EventArgs e)
        {
            var selNode = (tlvBranch)tlvGroupedFiles.GetModelObject(_clickedItem.Index);

            _numberNewNodes++;

            var ssg = new SpectrumSourceGroup();
            ssg.Name = (string.Format("{0}/New tlvBranch({1})", selNode.Text, _numberNewNodes).Replace(@"//",@"/"));

            var newNode = new tlvBranch
                              {
                                  Text = Path.GetFileName(ssg.Name),
                                  cms = cmRightClickGroupNode,
                                  Parent = selNode,
                                  Children = new List<object>(),
                                  Data = ssg
                              };

            selNode.Children.Add(newNode);
            OrganizeNode(selNode);
            tlvGroupedFiles.RefreshObject(selNode);
            if (selNode.Children.Count == 1)
                tlvGroupedFiles.Expand(selNode);
            tlvGroupedFiles.EditSubItem(tlvGroupedFiles.ModelToItem(newNode),0);
        }

        private void renameGroupToolStripMenuItem_Click(object sender, EventArgs e)
        {
            var selNode = (tlvBranch)tlvGroupedFiles.GetModelObject(_clickedItem.Index);


            if (selNode.Parent.Text != null)
                tlvGroupedFiles.EditSubItem(_clickedItem, 0);
        }

        private void removeGroupToolStripMenuItem_Click(object sender, EventArgs e)
        {
            var selNode = (tlvBranch)tlvGroupedFiles.GetModelObject(_clickedItem.Index);
            RemoveGroupNode(selNode);
        }

        private void RemoveGroupNode(tlvBranch selNode)
        {
            var abandonedSpectraSources = new List<tlvBranch>();
            getListOfSprectrumSourcesRecursively(selNode, ref abandonedSpectraSources);

            if (selNode.Parent.Text != null)
            {
                selNode.Parent.Children.Remove(selNode);
                tlvGroupedFiles.RefreshObject(selNode.Parent);
            }
            else
            {
                selNode.Children.Clear();
                tlvGroupedFiles.RefreshObject(selNode);
            }

            foreach (tlvBranch tn in abandonedSpectraSources)
            {
                var lvi = new ListViewItem {Text = Path.GetFileName(tn.Text), Tag = tn};
                lvNonGroupedFiles.Items.Add(lvi);
            }
        }

        private void RemoveFileNode(tlvBranch selNode)
        {
            selNode.Parent.Children.Remove(selNode);
            tlvGroupedFiles.RefreshObject(selNode.Parent);

            var lvi = new ListViewItem {Text = Path.GetFileName(selNode.Text), Tag = selNode};
            lvNonGroupedFiles.Items.Add(lvi);
        }

        private void removeFileNodeToolStripMenuItem_Click(object sender, EventArgs e)
        {
            var selNode = (tlvBranch)tlvGroupedFiles.GetModelObject(_clickedItem.Index);
            RemoveFileNode(selNode);
        }

        private void lvNonGroupedFiles_ItemDrag(object sender, ItemDragEventArgs e)
        {
            DoDragDrop(e.Item, DragDropEffects.Move);
        }

        private void lvNonGroupedFiles_DragDrop(object sender, DragEventArgs e)
        {
            if (e.Data.GetDataPresent("System.String", true))
            {
                var dragSources = from tlvBranch branch in tlvGroupedFiles.SelectedObjects
                                where !(branch.Data is SpectrumSourceGroup)
                                select branch;

                var remaining = from tlvBranch branch in tlvGroupedFiles.SelectedObjects
                                where (branch.Data is SpectrumSourceGroup)
                                select branch;
                
                foreach (var item in dragSources)
                    RemoveGroupNode(item);
                foreach (var item in remaining)
                    RemoveGroupNode(item);
                
                tlvGroupedFiles.Sort();
            }
        }

        private void lvNonGroupedFiles_DragEnter(object sender, DragEventArgs e)
        {
            e.Effect = DragDropEffects.Move;
        }

        private void lvNonGroupedFiles_DragOver(object sender, DragEventArgs e)
        {
            e.Effect = DragDropEffects.Move;
        }

        private void lvNonGroupedFiles_KeyDown(object sender, KeyEventArgs e)
        {
            if ((e.Modifiers & Keys.Control) == Keys.Control && e.KeyCode == Keys.A)
                foreach (ListViewItem lvi in lvNonGroupedFiles.Items)
                    lvi.Selected = true;
        }

        private void miResetFiles_Click(object sender, EventArgs e)
        {
            if (MessageBox.Show("Are you sure you want to remove all groups?","Remove All",MessageBoxButtons.YesNo) == DialogResult.Yes)
                RemoveGroupNode((tlvBranch)tlvGroupedFiles.GetModelObject(0));
        }

        private void miExpandGroups_Click(object sender, EventArgs e)
        {
            tlvGroupedFiles.ExpandAll();
        }

        private void miCollapseGroups_Click(object sender, EventArgs e)
        {
            tlvGroupedFiles.CollapseAll();
        }

        private void tlvGroupedFiles_CanDrop(object sender, OlvDropEventArgs e)
        {
            if (e.DataObject.GetType().ToString() == "System.Windows.Forms.DataObject")
                e.Effect = DragDropEffects.Move;
            else if (e.DataObject is OLVDataObject
                     && (e.DropTargetItem.RowObject is tlvBranch
                         && ((e.DropTargetItem.RowObject as tlvBranch).Data is SpectrumSourceGroup
                              || (e.DropTargetItem.RowObject as tlvBranch).Data is SpectrumSource)))
            {
                var target = e.DropTargetItem.RowObject as tlvBranch;
                var dragging = (e.DataObject as OLVDataObject).ModelObjects;
                var isValid = true;
                foreach (var item in dragging)
                {
                    if (checkIfDestGroupAChildNodeOfMe(target,item as tlvBranch))
                    {
                        isValid = false;
                        break;
                    }
                }
                if (isValid)
                    e.Effect = DragDropEffects.Move;
            }
        }

        private void tlvGroupedFiles_Dropped(object sender, OlvDropEventArgs e)
        {
            var target = (tlvBranch)e.DropTargetItem.RowObject;
            var index = -1;
            if (target.Data is SpectrumSource)
            {
                index = target.Parent.Children.IndexOf(target);
                target = target.Parent;
            }
            if (e.DataObject.GetType().ToString() == "System.Windows.Forms.DataObject")
            {
                var sources = from ListViewItem item in lvNonGroupedFiles.SelectedItems
                              select (tlvBranch) item.Tag;

                foreach (var source in sources)
                {
                    if (index >= 0)
                        target.Children.Insert(index, source);
                    else
                        target.Children.Add(source);
                    source.Parent = target;
                }

                tlvGroupedFiles.RefreshObject(target);

                var usedItems = from ListViewItem item in lvNonGroupedFiles.SelectedItems
                                select item;
                foreach (var item in usedItems)
                    lvNonGroupedFiles.Items.Remove(item);
            }
            else if (e.DataObject is OLVDataObject)
            {
                var dragging = e.DataObject as OLVDataObject;
                var sources = from tlvBranch item in dragging.ModelObjects
                              where item.Data is SpectrumSource
                              select item;
                var groups = from tlvBranch item in dragging.ModelObjects
                             where (item.Data is SpectrumSourceGroup
                             && (item.Data as SpectrumSourceGroup).Name != "\\")
                             select item;
                var sourcesToIgnore = new List<tlvBranch>();

                foreach (var group in groups)
                {
                    //find and ignore spectra in group
                    getListOfSprectrumSourcesRecursively(group, ref sourcesToIgnore);

                    group.Parent.Children.Remove(group);
                    tlvGroupedFiles.RefreshObject(group.Parent);
                    group.Parent = target;
                    if (target.Children.Any())
                        target.Children.Insert(0, group);
                    else
                        target.Children.Add(group);
                }

                sources = from tlvBranch s in sources where !sourcesToIgnore.Contains(s) select s;
                foreach (var source in sources)
                {
                    source.Parent.Children.Remove(source);
                    tlvGroupedFiles.RefreshObject(source.Parent);
                    source.Parent = target;
                    if (index >= 0)
                        target.Children.Insert(index, source);
                    else
                        target.Children.Add(source);
                }

                tlvGroupedFiles.RefreshObject(target);
                tlvGroupedFiles.Expand(target);
            }
            OrganizeNode(target);
        }

        private void tlvGroupedFiles_CellEditFinishing(object sender, CellEditEventArgs e)
        {
            if (e.NewValue.ToString().Contains("/"))
            {
                e.Cancel = true;
                return;
            }

            foreach (tlvBranch item in ((tlvBranch) e.RowObject).Parent.Children)
                if (item.Text == (string)e.NewValue)
                {
                    e.Cancel = true;
                    return;
                }
            ((tlvBranch)e.RowObject).Text = (string)e.NewValue;
        }

        private void tlvGroupedFiles_CellEditStarting(object sender, CellEditEventArgs e)
        {
            if (((tlvBranch)e.RowObject).Data is SpectrumSource ||(string)e.Value == "/")
            {
                e.Cancel = true;
            }

        }

        private void ApplyDefaultGroups(object sender, EventArgs e)
        {
            if (sender == miDefaultGroups &&
                MessageBox.Show("Are you sure you want to reset the groups to their initial values?",
                "Reset Groups?", MessageBoxButtons.YesNo) != DialogResult.Yes)
                return;

            var obj = tlvGroupedFiles.Objects;
            foreach(var item in obj)
                tlvGroupedFiles.RemoveObject(item);
            lvNonGroupedFiles.Items.Clear();

            var sourcesByGroup = session.Query<SpectrumSource>().Where(o => o.Group != null)
                                                                .ToLookup(o => (int)o.Group.Id.Value)
                                                                .ToDictionary(o => o.Key, o => o.ToList());

            var groups = session.Query<SpectrumSourceGroup>().ToList();
            var groupNode = new tlvBranch
            {
                Text = "/",
                Parent = new tlvBranch { Text = null },
                Children = new List<object>(),
                Data = (from g in groups where g.Name == "/" select g).Single()
            };

            groups.Remove(groupNode.Data as SpectrumSourceGroup);
            groupNode.Text = (groupNode.Data as SpectrumSourceGroup).Name;
            groupNode.cms = cmRightClickGroupNode;

            AddAssociatedSpectra(ref groupNode, sourcesByGroup);

            groupNode = FillBranch(groupNode, groups, sourcesByGroup);

            tlvGroupedFiles.AddObject(groupNode);
            _rootNode = groupNode;

            tlvGroupedFiles.ExpandAll();

            // add ungrouped sources
            foreach (var ss in session.Query<SpectrumSource>().Where(g => g.Group == null))
            {
                var lvi = new ListViewItem { Text = ss.Name, Tag = ss };
                lvNonGroupedFiles.Items.Add(lvi);
            }
        }

        private int OrganizeNode(tlvBranch target)
        {
            var sources = from tlvBranch item in target.Children
                           where item.Data is SpectrumSource
                           select item;
            var groups = from tlvBranch item in target.Children
                          where item.Data is SpectrumSourceGroup
                          select item;
            if (!sources.Any() && !groups.Any())
                return -1;
            
            target.Children = new List<object>();
            foreach (var item in groups)
                target.Children.Add(item);
            foreach (var item in sources)
                target.Children.Add(item);

            return groups.Count();
        }

    }
}