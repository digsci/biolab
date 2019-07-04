import pandas as pd

# Input files required

# All the 2016 data samples
data_file = "2016-data.csv"
# Lookup table from FeraID to CustomerID
lookuptable_file = "lookuptable.csv"
# Diatom abundances file from analysis pipeline
abundances_file = "Abundances.pass.csv"
# Output file name
darleq3_file = "darleq3-input.xlsx"

def seasontodate(season):
    if season == 'Spr':
        return '01/03/2016'
    elif season == 'Aut':
        return '01/09/2016'

def sample_data(feraID):
  sample = data_df[data_df['FERAID'] == feraID][['FERAID', 'BIOSYS site ID','Season','Alkalinity', 'SPT_Name']]
  return sample

def renameSampleIDs(abundances_df, lookup_df):
  def sampleID(customerID):
    sample_id = lookup_df[lookup_df['Customer_ID'] == customerID]['Sample_ID'].iat[0]
    return sample_id
  abundances_df.rename(mapper=sampleID, axis=1, inplace=True)
  return abundances_df

def getDataFrames():
  data_df = pd.read_csv(data_file)
  lookup_df = pd.read_csv(lookuptable_file)
  abundances_df = pd.read_csv(abundances_file, index_col=0, header=0, dtype=str)
  abundances_df = renameSampleIDs(abundances_df, lookup_df)
  return [data_df, lookup_df, abundances_df]


def createHeaders(abundances_df):
  header_df = pd.DataFrame()
  for sample_id in abundances_df:
    header_df = pd.concat([header_df,sample_data(int(sample_id))], ignore_index=True)
  header_df = header_df.T
  for sample in header_df:
    header_df[sample]['Season'] = seasontodate(header_df[sample]['Season'])
  header_df = header_df.rename({'FERAID':'SampleID', 'BIOSYS site ID':'SiteID', 'Season':'SampleDate','SPT_Name':'Reach'})
  return header_df

def createAbundance(abundance_df):
  abundance_df.columns = range(0, 200, 1)
  return abundance_df

def createDarleq3Input(header_df, abundance_df):
  header_df.insert(0,'','')
  abundance_df.insert(0,'',abundance_df.index)
  darleq3_input_df = pd.concat([header_df, abundance_df])
  return darleq3_input_df

if __name__ == "__main__":
  data_df, lookup_df, abundances_df = getDataFrames()
  header_df = createHeaders(abundances_df)
  abundance_df = createAbundance(abundances_df)
  darleq3_input_df = createDarleq3Input(header_df, abundance_df)
  darleq3_input_df.to_excel(darleq3_file)
  
