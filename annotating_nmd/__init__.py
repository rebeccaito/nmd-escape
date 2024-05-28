#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import warnings


class ParentWarnings(Warning):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


class IncorrectColumnWarning(ParentWarnings):
    pass


class HGVSpPatternWarning(ParentWarnings):
    pass


def check_bed(bed_df):
    """
    Checks whether dataframe is in 6-column bed format with expected headers. If the dataframe has fewer than 6 columns
    OR the dataframe has more than 6 columns, and the 6 expected column names are not present, then an Exception will be
    raised. If the dataframe is exactly 6 columns, but the column names are not as expected, a Warning is raised and the
    columns are renamed.

    :param bed_df: pandas dataframe
    :return: same dataframe or with renamed columns
    """
    col_names = ['chrom', 'start', 'end', 'cds_id', 'score', 'strand']

    if (len(bed_df.columns) < 6) or (len(bed_df.columns) > 6 and not all([col in bed_df.columns for col in col_names])):
        raise Exception("BED file is not in 6 column format")

    if list(bed_df.columns) != col_names and len(bed_df.columns) == 6:
        warnings.warn(f'Column names in bed dataframe do not match: {col_names}. Renaming now.', IncorrectColumnWarning)
        rename_dict = {}
        for i, key in enumerate(bed_df.columns):
            rename_dict[key] = col_names[i]
        bed_df.rename(columns=rename_dict, inplace=True)

    return bed_df


def sort_transcript_bed(transcript_df):
    """
    Helper function to check and sort transcript_df. Sorting by transcript strand is necessary for determining first
    vs. last coding exons.

    :param transcript_df: pandas dataframe with 1 transcript's coding exons in 6-column bed format
    :return: same dataframe sorted by gene orientation and by start coordinate
    """
    transcript_df = check_bed(transcript_df).copy()

    if 'cds_size' not in transcript_df.columns:
        transcript_df = preprocess_bed(transcript_df)

    if transcript_df.strand.values[0] == '+':
        ascending_sort = False
    else:
        ascending_sort = True

    transcript_df.sort_values('start', ascending=ascending_sort, inplace=True)
    transcript_df.reset_index(inplace=True)

    return transcript_df


def preprocess_bed(bed_df, capture_pattern='(^.*\.\d*)_cds.*$'):
    """
    Pre-process CDS bed dataframe to 1) extract transcript_name from CDS ID, and, 2) determine regions sizes.
    CDS ID defaults to be the format downloaded from UCSC Genome Browser for CDS bed for RefSeq transcripts

    :param bed_df: pandas dataframe with 6-columns bed format and the following named columns:
        ['chrom', 'start', 'end', 'cds_id', 'score', 'strand']
    :return: processed bed_df
    """

    bed_df = check_bed(bed_df)

    # drop dups due to PAR
    bed_df = bed_df.drop_duplicates(['cds_id', 'start'])

    bed_df['transcript_name'] = bed_df.cds_id.str.extract(capture_pattern)
    bed_df['cds_size'] = bed_df.end - bed_df.start

    return bed_df


def get_nmd_escape_size(transcript_df):
    """
    Determine NMD(-) size per transcript
    If transcript on + strand, grab the size of last (relative to genome) CDS + up to 55 bp of penultimate CDS
    If transcript on - strand, grab the size of first (relative to genome) CDS + up to 55 bp of second CDS

    :param transcript_df: CDS bed dataframe per transcript  with 6-columns bed format and an additional column for
            "cds_size" assuming preprocess_bed has been run
    :return: total size of NMD(-) region per transcript following 55 nt rule
    """

    transcript_df = sort_transcript_bed(transcript_df).copy()

    nmd_size = 0
    for i, row in transcript_df.iterrows():
        if i == 0:
            # if first exon in df (meaning last exon in transcript), grab entire CDS of exon
            nmd_size += row.cds_size
        elif i == 1:
            # if second exon in df (meaning penultimate exon in transcript), grab the minimum of 55nt or exon size
            nmd_size += min([row.cds_size, 55])
        else:
            return nmd_size

    return nmd_size


def get_nmd_escape_boundaries(transcript_df):
    """
    Create a bed dataframe of NMD(-) regions.
    If transcript is on minus strand, sort by ascending start coordinate, and grab first 2 CDS exons (last 2 in transcript)
    If transcript is on plus strand, sort by descending start coordinate, and grab first 2 CDS exons (last 2 in transcript)

    :param transcript_df: pandas dataframe with 6-columns bed format and the following named columns:
        ['chrom', 'start', 'end', 'cds_id', 'score', 'strand']
    :return: nmd_bed dataframe
    """

    transcript_df = sort_transcript_bed(transcript_df)

    nmd_bed = pd.DataFrame()
    for i, row in transcript_df.iterrows():
        if i == 0:
            nmd_bed = pd.concat([nmd_bed, pd.DataFrame(row).transpose()], ignore_index = True)
        elif i == 1:
            remaining_nmd_size = min([55, row.cds_size]) # if cds is bigger than needed, then adjust start/stop to be size 55
            if row.strand == '+':
                row['start'] = row.end - remaining_nmd_size
            else:
                row['end'] = row.start + remaining_nmd_size
            row['cds_size'] = row.end - row.start
            nmd_bed = pd.concat([nmd_bed, pd.DataFrame(row).transpose()], ignore_index = True)
        else:
            return nmd_bed

    return nmd_bed


def make_boundaries_df(bed_df):
    """
    Convenience wrapper to apply get_nmd_escape_boundaries to ALL transcripts in a pandas dataframe

    :param bed_df: pandas dataframe with 6-columns bed format and the following named columns:
        ['chrom', 'start', 'end', 'cds_id', 'score', 'strand']
    :return: boundaries_df, pandas dataframe in bed format with NMD(-) boundaries
    """
    if 'transcript_name' not in bed_df.columns:
        bed_df = preprocess_bed(bed_df)

    boundaries_df = bed_df.groupby('transcript_name').apply(get_nmd_escape_boundaries, include_groups=False)
    return boundaries_df


def make_cds_size_df(bed_df):
    """
    Convenience function to take a CDS bed dataframe and return total CDS length, NMD(-) length, PDOT length, and
        pdot start position for the NMD(-) region.

    :param bed_df: pandas dataframe with 6-columns bed format and the following named columns:
        ['chrom', 'start', 'end', 'cds_id', 'score', 'strand']
    :return: sizes: pandas dataframe with per transcript CDS, NMD(-), and pdot lengths
    """

    bed_df = check_bed(bed_df)

    if 'cds_size' not in bed_df.columns:
        bed_df = preprocess_bed(bed_df)

    cds_sizes = bed_df.groupby('transcript_name').cds_size.sum().reset_index().rename(columns={0: 'cds_size'})
    nmd_sizes = bed_df.groupby('transcript_name').apply(
        get_nmd_escape_size, include_groups=False
    ).reset_index().rename(columns={0: 'nmd_escape_size'})

    # combine cds size and nmd(-) size into 1 dataframe
    sizes = cds_sizes.merge(nmd_sizes, on='transcript_name', how='outer')

    # calculate pdot for NMD(-)
    sizes['total_pdot_length'] = (sizes.cds_size / 3).astype(int)
    sizes['nmd_pdot_start'] = sizes.total_pdot_length - (sizes.nmd_escape_size / 3)

    return sizes


def get_upstream_frameshift(annotated_df, nmd_df):
    """
        This function will return whether a frameshift variant is PTVesc
        even if the variant is upstream of NMD(-) as long as the stop codon
        falls within NMD(-)

    :param annotated_df: pandas dataframe with variant annotations including "HGVSp".
    :param nmd_df:
    :return:
    """

    hgvsp_pattern = '^.*\.\d*:p\..{3}(\d*).{3}fsTer(\d*)'

    annotated_df[['var_pdot', 'stop_pdot_shift']] = annotated_df.HGVSp.str.extract(hgvsp_pattern)
    annotated_df['stop_pdot_shift'] = annotated_df.stop_pdot_shift.replace("^$", np.nan, regex=True)

    # add a check for matching expected HGVSp pattern
    if (
            (annotated_df.var_pdot.isna())
            | (annotated_df.stop_pdot_shift.isna())
    ).all():
        warnings.warn("No variants meet the expected pattern for HGVSp truncating frameshift variants.\n"
                "Please refer to https://varnomen.hgvs.org/recommendations/protein/variant/frameshift/",
                      HGVSpPatternWarning)

    annotated_df['stop_pdot'] = annotated_df.var_pdot.astype('Int64') + annotated_df.stop_pdot_shift.astype('Int64')
    annotated_df = annotated_df.merge(nmd_df[['transcript_name', 'nmd_pdot_start']],
                                   on='transcript_name', how='left')
    annotated_df['is_nmd_frameshift'] = annotated_df.stop_pdot > annotated_df.nmd_pdot_start

    return annotated_df
