export default {
  multipass: true,
  plugins: [
    {
      name: 'preset-default',
      params: {
        overrides: {
          collapseGroups: false,
          mergePaths: false
        }
      }
    }
  ]
};
