
import pandas as pd
import numpy as np
import seaborn as sns
from collections import OrderedDict
from functools import reduce
import json
from scipy import interpolate
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

class Demog:

    name_conversion_dict = {'Age (x)': 'Age',
                            'Central death rate m(x,n)': 'Mortality_mid',
                            'Age interval (n)': 'Interval',
                            'Period': 'Years'
                            }

    sex_dict = {'Male': 0, 'Female': 1}

    @staticmethod
    def min_not_nan(x_list):
        loc_in = list(filter(lambda x: not np.isnan(x), x_list))
        return np.min(loc_in)

    @staticmethod
    def map_year(x_tuple, flag='mid'):

        valid_entries_loc = ['mid', 'end', 'start']

        if flag not in valid_entries_loc:
            raise ValueError('invalid endpoint specified')

        if flag == 'mid':
            return (x_tuple[0] + x_tuple[1]) / 2.0
        elif flag == 'start':
            return x_tuple[0]
        else:
            return x_tuple[1]

    @staticmethod
    def construct_interval(x, y):
        return (x, x + y)

    @staticmethod
    def midpoint(x, y):
        return (x + y) / 2.0

    @staticmethod
    def generate_dict_order(tuple_list, which_entry=1):

        my_unordered_list = tuple_list.apply(lambda x: x[which_entry])
        dict_to_order = OrderedDict(zip(tuple_list, my_unordered_list))

        return dict_to_order

    def __init__(self, file_male,file_female, regression_interval=[1970, 1980],which_point = 'mid', predict_horizon=2062.5):
        self.filenames = [file_male, file_female]
        self.interval_fit = regression_interval
        self.predict_horizon = predict_horizon
        self.estimates_list = []
        self.which_point = which_point

    def do_extrapolation(self, add_to_list=True):
        df_mort_male = pd.read_csv(self.filenames[0], usecols=Demog.name_conversion_dict)
        df_mort_male['Sex'] = 'Male'
        df_mort_female = pd.read_csv(self.filenames[1], usecols=Demog.name_conversion_dict)
        df_mort_female['Sex'] = 'Female'
        df_mort = pd.concat([df_mort_male, df_mort_female], axis=0)
        df_mort.rename(columns=Demog.name_conversion_dict, inplace=True)
        df_mort['Years'] = df_mort['Years'].apply(lambda xx: tuple(
            [float(zz) for zz in xx.split('-')]))  # this might be a bit too format specific (ie dashes in input)

        # log transform the data and drop unneeded columns
        df_mort['log_Mortality_mid'] = df_mort['Mortality_mid'].apply(lambda x: np.log(x))
        df_mort['Age'] = df_mort[['Age', 'Interval']].apply(lambda zz: Demog.construct_interval(*zz), axis=1)

        year_order_dict = Demog.generate_dict_order(df_mort['Years'])
        age_order_dict = Demog.generate_dict_order(df_mort['Age'])
        df_mort['sortby2'] = df_mort['Age'].map(age_order_dict)
        df_mort['sortby1'] = df_mort['Sex'].map(Demog.sex_dict)
        df_mort['sortby3'] = df_mort['Years'].map(year_order_dict)
        df_mort.sort_values(['sortby1', 'sortby2', 'sortby3'], inplace=True)
        df_mort.drop(columns=['Mortality_mid', 'Interval', 'sortby1', 'sortby2', 'sortby3'], inplace=True)

        # convert to years (and to string for age_list due to really annoying practical slicing reasons
        df_mort['Years'] = df_mort['Years'].apply(lambda x: Demog.map_year(x, self.which_point))
        df_mort['Age'] = df_mort['Age'].apply(lambda x: str(x))
        df_before_time = df_mort[df_mort['Years'].between(0, self.interval_fit[0])].copy()

        df_mort.set_index(['Sex', 'Age'], inplace=True)
        sex_list = list(set(df_mort.index.get_level_values('Sex')))
        age_list = list(set(df_mort.index.get_level_values('Age')))

        df_list = []
        df_list_future = []
        for sex in sex_list:
            for age in age_list:
                tmp_data = df_mort.loc[(sex, age, slice(None)), :]
                extrap_model = LinearRegression(normalize=True, fit_intercept=True, copy_X=True)
                extrap_model_2 = LinearRegression(normalize=True, fit_intercept=True, copy_X=True)

                first_extrap_df = tmp_data[tmp_data['Years'].between(self.interval_fit[0], self.interval_fit[1])]
                XX = tmp_data[tmp_data['Years'].between(self.interval_fit[0], self.predict_horizon)].values[:, 0]

                values = first_extrap_df.values
                extrap_model.fit(values[:, 0].reshape(-1, 1), values[:, 1])

                extrap_predictions = extrap_model.predict(XX.reshape(-1, 1))

                loc_df = pd.DataFrame.from_dict({'Sex': sex, 'Age': age, 'Years': XX, 'Extrap': extrap_predictions})
                loc_df.set_index(['Sex', 'Age', 'Years'], inplace=True)

                df_list.append(loc_df.copy())

        df_e1 = pd.concat(df_list, axis=0)

        df_list_final = [df_mort, df_e1]
        df_total = reduce(lambda left, right: pd.merge(left, right, on=['Sex', 'Age', 'Years']), df_list_final)

        df_total = df_total.reset_index(inplace=False).set_index(['Sex', 'Age'], inplace=False)

        df_total['Extrap'] = df_total['Extrap'].apply(np.exp)
        df_total['Data'] = df_total['log_Mortality_mid'].apply(np.exp)
        df_before_time['Data'] = df_before_time['log_Mortality_mid'].apply(np.exp)

        df_before_time.set_index(['Sex', 'Age'], inplace=True)
        df_total = pd.concat([df_total, df_before_time], axis=0, join='outer', sort=True)

        df_total.reset_index(inplace=True)
        df_total['sortby2'] = df_total['Age'].map(age_order_dict)
        df_total['sortby1'] = df_total['Sex'].map(Demog.sex_dict)
        df_total.sort_values(by=['sortby1', 'sortby2', 'Years'], inplace=True)
        df_total.drop(columns=['sortby1', 'sortby2'], inplace=True)

        if not add_to_list:
            self.estimates_list = [df_total.copy()]
        else:
            self.estimates_list.append(df_total.copy())

    def makeplots(self, save=None):

        for i,df_loc in enumerate(self.estimates_list):

            df_loc.reset_index(inplace=True)

            grid = sns.FacetGrid(df_loc, col='Age', row='Sex')
            grid.map(plt.plot, 'Years', 'Extrap', color='b', label='Pre-epidemic')
            grid.map(plt.plot, 'Years', 'Data', color='r', linestyle='--', label='Data', marker='None')
            df_loc['Min_Mortality'] = df_loc[['Data', 'Extrap']].apply(Demog.min_not_nan, axis=1)

            if save:
                plt.savefig(str(save) + '_' + str(i) + '.png')
            else:
                plt.show()

    def create_json_overlay(self, template, output_name='Extract_demog.json', csv_out=False, n=0, results_scale_factor=0.001/365.0):

        df = self.estimates_list[n]
        df['FE'] = df[['Data', 'Extrap']].apply(Demog.min_not_nan, axis=1)
        df['Age'] = df['Age'].apply(lambda x: int(x.split(',')[1].split(')')[0]))
        male_df = df[df['Sex'] == 'Male']
        female_df = df[df['Sex'] == 'Female']

        male_df.set_index(['Sex','Age', 'Years'], inplace=True)
        female_df.set_index(['Sex','Age','Years'], inplace=True)
        male_data = male_df['FE']
        female_data = female_df['FE']

        male_data = male_data.unstack(-1)
        male_data.sort_index(level ='Age', inplace=True)
        female_data = female_data.unstack(-1)
        female_data.sort_index(level='Age', inplace=True)

        years_out_male = list(male_data.columns)
        years_out_female =list( female_data.columns)

        age_out_male = list(male_data.index.get_level_values('Age'))
        age_out_female = list(male_data.index.get_level_values('Age'))

        male_output= male_data.values
        female_output = female_data.values

        if csv_out:
            male_data.to_csv('Male' + csv_out)
            female_data.to_csv('Female' + csv_out)

        with open(template,'r') as f:
            out_json = json.load(f)

        dict_female = {'NumPopulationGroups': list(female_data.shape),
                       'AxisNames': ['age', 'year'],
                       'AxisScaleFactors': [365.0, 1],
                       'AxisUnits': ['years', 'years'],
                       'NumDistributionAxes': 2,
                       'PopulationGroups': [age_out_female, years_out_female],
                       'ResultsScaleFactor': results_scale_factor,
                       'ResultsUnits': 'annual deaths per capita',
                       'ResultsValues': female_output.tolist()
                       }

        dict_male =   {'NumPopulationGroups': list(male_data.shape),
                       'AxisNames': ['age', 'year'],
                       'AxisScaleFactors': [365.0, 1],
                       'AxisUnits': ['years', 'years'],
                       'NumDistributionAxes': 2,
                       'PopulationGroups': [age_out_male, years_out_male],
                       'ResultsScaleFactor': results_scale_factor,
                       'ResultsUnits': 'annual deaths per capita',
                       'ResultsValues': male_output.tolist()
                       }
        out_json['Defaults']['IndividualAttributes']['MortalityDistributionFemale'] = dict_female
        out_json['Defaults']['IndividualAttributes']['MortalityDistributionMale'] = dict_male

        out_data = json.dumps(out_json, indent=4)

        with open(output_name,'w') as f_out:
            f_out.write(out_data)




if __name__ == '__main__':
    # placeholders for command line arguments
    make_plots = True

    # User variables where to start extrapolation (ie before HIV epidemic) and how far to look back
    # All in years unless noted otherwise
    extrap_point = 1980
    look_back = 20
    time_horizon = 2100
    now_ish = 1980
    look_back_from_now_ish = 20

    filename_male = 'Malawi_male_mortality.csv'
    filename_female = 'Malawi_female_mortality.csv'

    interval_for_fitting_regression = [1970, 1980]
    predict_horizons = 2062.5
    #predict_horizons = 2050

    myob = Demog(filename_male, filename_female) #creatrs object
    myob.do_extrapolation()  #does extrapolation (if you do more than one it adds it to a list so you can compare)
    myob.makeplots(save=False) # make plots by age
    myob.create_json_overlay('template.json',output_name='Output_demog.json', csv_out='Test.csv') #create json
