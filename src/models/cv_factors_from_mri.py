import logging

import click
import click_log
import numpy as np
import sklearn.dummy
import sklearn.ensemble
import sklearn.metrics
import sklearn.model_selection
import sklearn.neighbors
import sklearn.pipeline
import sklearn.preprocessing
import xarray as xr


logger = logging.getLogger(__name__)


models = {
    'mean': sklearn.pipeline.Pipeline([
        ('impute', sklearn.preprocessing.Imputer(strategy='median', axis=0)),
        ('dummy', sklearn.dummy.DummyRegressor('mean')),
    ]),
    'random_forest': sklearn.pipeline.Pipeline([
        ('impute', sklearn.preprocessing.Imputer()),
        ('rf', sklearn.ensemble.RandomForestRegressor(
            n_estimators=1000,
            n_jobs=4,
            max_depth=2,
        )),
    ]),
    'knn': sklearn.pipeline.Pipeline([
        ('impute', sklearn.preprocessing.Imputer()),
        ('scale', sklearn.preprocessing.RobustScaler()),
        ('select_features', sklearn.feature_selection.SelectKBest(
            sklearn.feature_selection.f_regression, k=4,
        )),
        ('knn', sklearn.neighbors.KNeighborsRegressor(
            n_neighbors=20,
            n_jobs=4,
        )),
    ]),
}


def named_full(coords, *args, **kwargs):
    shape = [len(c[1]) for c in coords]
    return xr.DataArray(
        np.full(shape, *args, **kwargs),
        dims=[c[0] for c in coords],
        coords={c[0]: c[1] for c in coords},
    )


def cv(model, X, y):
    cv = sklearn.model_selection.KFold(10)
    y_predicted = np.full(y.shape, np.nan, y.dtype)
    for train_idx, test_idx in cv.split(X):
        model.fit(X[train_idx, :], y[train_idx])
        y_predicted[test_idx] = model.predict(X[test_idx, :])
    return y_predicted


def cv_all_models(mri, sfa):
    logger.info("Cross Validating")
    mse_cv = named_full([
        ('model', np.array(list(models.keys()), object)),
        ('factor', sfa.coords['factor']),
    ], np.nan)
    ev_cv = named_full([
        ('model', np.array(list(models.keys()), object)),
        ('factor', sfa.coords['factor']),
    ], np.nan)
    mae_cv = named_full([
        ('model', np.array(list(models.keys()), object)),
        ('factor', sfa.coords['factor']),
    ], np.nan)
    feature_importance = named_full([
        ('model', np.array(list(models.keys()), object)),
        ('factor', sfa.coords['factor']),
        ('cad_feature', mri.coords['cad_feature']),
    ], np.nan)

    for factor in sfa['factor']:
        for model_id, model in models.items():
            logger.info("\tCV {} {}".format(model_id, str(factor.values)))
            X = mri.values
            y = sfa.loc[:, factor].values
            y_pred = cv(model, X, y)
            mse_cv.loc[model_id, factor] =\
                sklearn.metrics.mean_squared_error(y, y_pred)
            mae_cv.loc[model_id, factor] =\
                sklearn.metrics.median_absolute_error(y, y_pred)
            ev_cv.loc[model_id, factor] =\
                sklearn.metrics.explained_variance_score(y, y_pred)

            model.fit(X, y)
            if hasattr(model, 'feature_importances_'):
                feature_importance.loc[model_id, factor, :] =\
                    model.feature_importances_
            elif hasattr(model, 'steps'):
                if hasattr(model.steps[-1][1], 'feature_importances_'):
                    feature_importance.loc[model_id, factor, :] =\
                        model.steps[-1][1].feature_importances_

    return xr.Dataset({
        'mean_square_error': mse_cv,
        'median_absolute_error': mae_cv,
        'explained_variance': ev_cv,
        'feature_importance': feature_importance,
    })


@click.command()
@click.argument('mri', type=click.Path(exists=True, dir_okay=False,
                                       resolve_path=True))
@click.argument('sfa', type=click.Path(exists=True, dir_okay=False,
                                       resolve_path=True))
@click.argument('out', type=click.Path(exists=False, dir_okay=False,
                                       writable=True, resolve_path=True))
@click_log.simple_verbosity_option()
@click_log.init(__name__)
def cli(mri, sfa, out):
    """Train machine learning models to predict MRI features from factors."""
    logger.info("Loading SFA model")
    sfa_ds = xr.open_dataset(sfa).load()
    sfa_ds = sfa_ds.transpose('case', 'factor')
    sfa = sfa_ds['factors'].transpose('case', 'factor')

    logger.info("Loading MRI features")
    mri_ds = xr.open_dataset(mri).load()
    mri_ds = mri_ds.reindex(case=sfa['case'])
    mri_features = list(set(mri_ds.keys()) - {'case', 'Comment', 'MultiFocal'})
    mri = mri_ds[mri_features].to_array('cad_feature')
    mri = mri.transpose('case', 'cad_feature')

    if not all(mri['case'].values == sfa['case'].values):
        logger.error("MRI and SFA data do not match up.")
        exit(-1)

    cv_ds = cv_all_models(mri, sfa)
    cv_ds['factor_name'] = sfa_ds['factor_name']
    cv_ds.to_netcdf(out)


if __name__ == '__main__':
    cli()
